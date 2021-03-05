module SSAmod

export SSA, SSA_deg_zero, SSA_deg_zero_log, SSA_deg_zero_all_log, SSA_deg_zero_conc

using Distributions, LambertW;


function SSA(S_time::Int64, inf_params::Array{Float64,1}, spec_params::Array{Float64,1}, max_time::Int64)
""" function with 5 params to be inferred {a0,a1,b0,b1,d} """
    M = 2; N = 1; 
    α0 = inf_params[1]; α1 = inf_params[2]; β0 = inf_params[3]; β1 = inf_params[4]; d = inf_params[5]; 
    V0 = spec_params[1]; n0 = spec_params[2]; θ = spec_params[3];
    
    # Define stoichiometric matrix
    S_mat = zeros(M);
    S_mat[1] = +1; # the burst size is randomly chosen later.
    S_mat[2] = -1;
    
    sp = 1.0; # sampling period is same as experimental data.
    n = zeros(S_time,Int64(floor(max_time/sp))); # reaction state vector
    c = zeros(S_time,Int64(floor(max_time/sp))); # conc vector
    v = ones(S_time,Int64(floor(max_time/sp))); # volume vector
    
    """
    * This first loop is for if we decide to average over more than one lineage - sim = 1 by default.
    * It also defines some constants that will be needed for the sim.
    """
    for sim in 1:S_time
        n_temp = n0; # set the current state vec to appropriate molecule number at homeostasis.
        T = 0.0 ::Float64; # set the overall time

        # while loop choosing to stop the simulation once the maximum time is reached.
        while T <= max_time
            # incorporate the time dependent rate explicitly.

            # Step 1: Calculate tau chosen rection for time dep rate.
            r1,r2 = rand(1,2);
            a1 = α0*(V0*exp(θ*T))^α1; a2 = d*n_temp;
            if isinf(lambertw((a1/a2)*(r1^(-α1*θ/a2))*exp(a1/a2))) # (**) if value is infinite go again.
                # print(a1/a2, "\n", r1^(-α1/a2), "\n", exp(a1/a2),"\n")
                tau = 0.0; # set to zero by default if get error.
            elseif a2 > 0 # distinguish between a2=0 and a2>0 cases separately.
                # print("hi")
                tau = (1/(α1*θ*a2))*(a1+α1*θ*log(1/r1)-a2*lambertw((a1/a2)*(r1^(-α1*θ/a2))*exp(a1/a2)));
            else
                tau = (1/(α1*θ))*log(1+(α1*θ/a1)*log(1/r1));
            end
            tau_props = [(a1/(α1*θ))*(exp(α1*θ*tau)-1),tau*a2]; sum_props = sum(tau_props);
            next_r = findfirst(x -> x>=r2*sum_props,cumsum(tau_props));

            # Step 2: Update the system
            T += tau;
            if isnan((T-tau)/sp) # catch the numerical error from (**)
                print("error!", "\n", r1, "\n", α1/a2, "\n", r1^(-α1/a2), "\n", exp(a1/a2), "\n", exp(a1/a2), "\n", (a1/a2)*(r1^(-α1/a2))*exp(a1/a2), "\n","boom")
            end
            # Update the trajectory vector
            if T < max_time # checks that T hasen't exceeded cell-cycle time or total time.
                for t in Int(ceil((T-tau)/sp)) : Int(floor(T/sp))
                    n[sim,t+1] = n_temp;
                    c[sim,t+1] = n_temp/(V0*exp(θ*t)); # the (t*sp)-(T-tc) updates is the volume for that time indexed by t.
                    v[sim,t+1] = V0*exp(θ*t);
                end
                # Fire reaction next_r
                if next_r == 1
                    b = β0*(V0*exp(θ*T))^β1;
                    bs = rand(Geometric((1/b)/((1/b)+1))); # rand burst size.
                    prod = bs*S_mat[next_r];
                elseif next_r == 2
                    prod = S_mat[next_r];
                else
                    println("error!")
                end
                n_temp = n_temp + prod;
            else # else the cell cycle time alone has been exceeded and do this
                for t in Int(ceil((T-tau)/sp)) : Int(floor(max_time/sp)-1)
                    n[sim,t+1] = n_temp;
                    c[sim,t+1] = n_temp/(V0*exp(θ*t));
                    v[sim,t+1] = V0*exp(θ*t);
                end
            end
        end
        #if mod(sim,100) == 0 # for multiple sims keep output after every 100.
        #    println(sim)
        end
    return n,c,v
end


function SSA_deg_zero(S_time::Int64, inf_params::Array{Float64,1}, spec_params::Array{Float64,1}, max_time::Int64)
""" function with 4 params to be inferred {a0,a1,b0,b1}, assume d = 0. """
    M = 1; N = 1; 
    α0 = inf_params[1]; α1 = inf_params[2]; β0 = inf_params[3]; β1 = inf_params[4];
    V0 = spec_params[1]; n0 = spec_params[2]; θ = spec_params[3];
    
    # Define stoichiometric matrix
    S_mat = zeros(M);
    S_mat[1] = +1; # the burst size is randomly chosen later.
    
    sp = 1.0; # sampling period is same as experimental data.
    n = zeros(S_time,Int64(floor(max_time/sp))); # reaction state vector
    c = zeros(S_time,Int64(floor(max_time/sp))); # conc vector
    v = ones(S_time,Int64(floor(max_time/sp))); # volume vector
    
    """
    * This first loop is for if we decide to average over more than one lineage - sim = 1 by default.
    * It also defines some constants that will be needed for the sim.
    """
    for sim in 1:S_time
        n_temp = n0; # set the current state vec to appropriate molecule number at homeostasis.
        T = 0.0 ::Float64; # set the overall time

        # while loop choosing to stop the simulation once the maximum time is reached.
        while T <= max_time
            # incorporate the time dependent rate explicitly.

            # Step 1: Calculate tau chosen rection for time dep rate.
            r1 = rand(1)[1];
            a1 = α0*(V0*exp(θ*T))^α1;
            tau = (1/(α1*θ))*log(1+(α1*θ/a1)*log(1/r1));

            # Step 2: Update the system
            T += tau;
            # Update the trajectory vector
            if T < max_time # checks that T hasen't exceeded cell-cycle time or total time.
                for t in Int(ceil((T-tau)/sp)) : Int(floor(T/sp))
                    n[sim,t+1] = n_temp;
                    c[sim,t+1] = n_temp/(V0*exp(θ*t)); # the (t*sp)-(T-tc) updates is the volume for that time indexed by t.
                    v[sim,t+1] = V0*exp(θ*t);
                end
                # Fire reaction
                b = β0*(V0*exp(θ*T))^β1;
                bs = rand(Geometric((1/b)/((1/b)+1))); # rand burst size.
                prod = bs*S_mat[1];
                n_temp = n_temp + prod;
            else # else the cell cycle time alone has been exceeded and do this
                for t in Int(ceil((T-tau)/sp)) : Int(floor(max_time/sp)-1)
                    n[sim,t+1] = n_temp;
                    c[sim,t+1] = n_temp/(V0*exp(θ*t));
                    v[sim,t+1] = V0*exp(θ*t);
                end
            end
        end
        #if mod(sim,100) == 0 # for multiple sims keep output after every 100.
        #    println(sim)
        end
    return n,c,v
end

function SSA_deg_zero_log(S_time::Int64, inf_params::Array{Float64,1}, spec_params::Array{Float64,1}, max_time::Int64)
""" function with 4 params to be inferred {a0,a1,b0,b1}, assume d = 0. 
    some scaling params to be inferred are done so in log scale."""
    M = 1; N = 1; 
    α0 = 10^inf_params[1]; α1 = inf_params[2]; β0 = 10^inf_params[3]; β1 = inf_params[4];
    V0 = spec_params[1]; n0 = spec_params[2]; θ = spec_params[3];
    
    # Define stoichiometric matrix
    S_mat = zeros(M);
    S_mat[1] = +1; # the burst size is randomly chosen later.
    
    sp = 1.0; # sampling period is same as experimental data.
    n = zeros(S_time,Int64(floor(max_time/sp))); # reaction state vector
    c = zeros(S_time,Int64(floor(max_time/sp))); # conc vector
    v = ones(S_time,Int64(floor(max_time/sp))); # volume vector
    
    """
    * This first loop is for if we decide to average over more than one lineage - sim = 1 by default.
    * It also defines some constants that will be needed for the sim.
    """
    for sim in 1:S_time
        n_temp = n0; # set the current state vec to appropriate molecule number at homeostasis.
        T = 0.0 ::Float64; # set the overall time

        # while loop choosing to stop the simulation once the maximum time is reached.
        while T <= max_time
            # incorporate the time dependent rate explicitly.

            # Step 1: Calculate tau chosen rection for time dep rate.
            r1 = rand(1)[1];
            a1 = α0*(V0*exp(θ*T))^α1;
            tau = (1/(α1*θ))*log(1+(α1*θ/a1)*log(1/r1));

            # Step 2: Update the system
            T += tau;
            # Update the trajectory vector
            if T < max_time # checks that T hasen't exceeded cell-cycle time or total time.
                for t in Int(ceil((T-tau)/sp)) : Int(floor(T/sp))
                    n[sim,t+1] = n_temp;
                    c[sim,t+1] = n_temp/(V0*exp(θ*t)); # the (t*sp)-(T-tc) updates is the volume for that time indexed by t.
                    v[sim,t+1] = V0*exp(θ*t);
                end
                # Fire reaction
                b = β0*(V0*exp(θ*T))^β1;
                bs = rand(Geometric((1/b)/((1/b)+1))); # rand burst size.
                prod = bs*S_mat[1];
                n_temp = n_temp + prod;
            else # else the cell cycle time alone has been exceeded and do this
                for t in Int(ceil((T-tau)/sp)) : Int(floor(max_time/sp)-1)
                    n[sim,t+1] = n_temp;
                    c[sim,t+1] = n_temp/(V0*exp(θ*t));
                    v[sim,t+1] = V0*exp(θ*t);
                end
            end
        end
        #if mod(sim,100) == 0 # for multiple sims keep output after every 100.
        #    println(sim)
        end
    return n,c,v
end


function SSA_deg_zero_all_log(S_time::Int64, inf_params::Array{Float64,1}, spec_params::Array{Float64,1}, max_time::Int64)
""" function with 4 params to be inferred {a0,a1,b0,b1}, assume d = 0. 
    ALL scaling params to be inferred are done so in log scale."""
    M = 1; N = 1; 
    α0 = exp(inf_params[1]); α1 = exp(inf_params[2]); β0 = exp(inf_params[3]); β1 = exp(inf_params[4]);
    V0 = spec_params[1]; n0 = spec_params[2]; θ = spec_params[3];
    
    # Define stoichiometric matrix
    S_mat = zeros(M);
    S_mat[1] = +1; # the burst size is randomly chosen later.
    
    sp = 1.0; # sampling period is same as experimental data.
    n = zeros(S_time,Int64(floor(max_time/sp))); # reaction state vector
    c = zeros(S_time,Int64(floor(max_time/sp))); # conc vector
    v = ones(S_time,Int64(floor(max_time/sp))); # volume vector
    
    """
    * This first loop is for if we decide to average over more than one lineage - sim = 1 by default.
    * It also defines some constants that will be needed for the sim.
    """
    for sim in 1:S_time
        n_temp = n0; # set the current state vec to appropriate molecule number at homeostasis.
        T = 0.0 ::Float64; # set the overall time

        # while loop choosing to stop the simulation once the maximum time is reached.
        while T <= max_time
            # incorporate the time dependent rate explicitly.

            # Step 1: Calculate tau chosen rection for time dep rate.
            r1 = rand(1)[1];
            a1 = α0*(V0*exp(θ*T))^α1;
            tau = (1/(α1*θ))*log(1+(α1*θ/a1)*log(1/r1));

            # Step 2: Update the system
            T += tau;
            # Update the trajectory vector
            if T < max_time # checks that T hasen't exceeded cell-cycle time or total time.
                for t in Int(ceil((T-tau)/sp)) : Int(floor(T/sp))
                    n[sim,t+1] = n_temp;
                    c[sim,t+1] = n_temp/(V0*exp(θ*t)); # the (t*sp)-(T-tc) updates is the volume for that time indexed by t.
                    v[sim,t+1] = V0*exp(θ*t);
                end
                # Fire reaction
                b = β0*(V0*exp(θ*T))^β1;
                bs = rand(Geometric((1/b)/((1/b)+1))); # rand burst size.
                prod = bs*S_mat[1];
                n_temp = n_temp + prod;
            else # else the cell cycle time alone has been exceeded and do this
                for t in Int(ceil((T-tau)/sp)) : Int(floor(max_time/sp)-1)
                    n[sim,t+1] = n_temp;
                    c[sim,t+1] = n_temp/(V0*exp(θ*t));
                    v[sim,t+1] = V0*exp(θ*t);
                end
            end
        end
        #if mod(sim,100) == 0 # for multiple sims keep output after every 100.
        #    println(sim)
        end
    return n,c,v
end

""" Below here are all the inference method where the production and burst size are instead dependent on the protein concentrations"""
function SSA_deg_zero_conc(S_time::Int64, inf_params::Array{Float64,1}, spec_params::Array{Float64,1}, max_time::Int64)
""" function with 4 params to be inferred {a0,a1,b0,b1}, assume d = 0. """
    M = 1; N = 1; 
    α0 = inf_params[1]; α1 = inf_params[2]; β0 = inf_params[3]; β1 = inf_params[4];
    V0 = spec_params[1]; n0 = spec_params[2]; θ = spec_params[3];
    
    # Define stoichiometric matrix
    S_mat = zeros(M);
    S_mat[1] = +1; # the burst size is randomly chosen later.
    
    sp = 1.0; # sampling period is same as experimental data.
    n = zeros(S_time,Int64(floor(max_time/sp))); # reaction state vector
    c = zeros(S_time,Int64(floor(max_time/sp))); # conc vector
    v = ones(S_time,Int64(floor(max_time/sp))); # volume vector
    
    """
    * This first loop is for if we decide to average over more than one lineage - sim = 1 by default.
    * It also defines some constants that will be needed for the sim.
    """
    for sim in 1:S_time
        n_temp = n0; # set the current state vec to appropriate molecule number at homeostasis.
        T = 0.0 ::Float64; # set the overall time

        # while loop choosing to stop the simulation once the maximum time is reached.
        while T <= max_time
            # incorporate the time dependent rate explicitly.

            # Step 1: Calculate tau chosen rection for time dep rate.
            r1 = rand(1)[1];
            a1 = α0 * n_temp^α1 * (V0*exp(θ*T))^(-α1);
            tau = (1/(α1*θ))*log(a1/(a1+α1*θ*log(r1)));

            # Step 2: Update the system
            T += tau;
            # Update the trajectory vector
            if T < max_time # checks that T hasen't exceeded cell-cycle time or total time.
                for t in Int(ceil((T-tau)/sp)) : Int(floor(T/sp))
                    n[sim,t+1] = n_temp;
                    c[sim,t+1] = n_temp/(V0*exp(θ*t)); # the (t*sp)-(T-tc) updates is the volume for that time indexed by t.
                    v[sim,t+1] = V0*exp(θ*t);
                end
                # Fire reaction
                b = β0* n_temp^β1 *(V0*exp(θ*T))^(-β1);
                bs = rand(Geometric((1/b)/((1/b)+1))); # rand burst size.
                prod = bs*S_mat[1];
                n_temp = n_temp + prod;
            else # else the cell cycle time alone has been exceeded and do this
                for t in Int(ceil((T-tau)/sp)) : Int(floor(max_time/sp)-1)
                    n[sim,t+1] = n_temp;
                    c[sim,t+1] = n_temp/(V0*exp(θ*t));
                    v[sim,t+1] = V0*exp(θ*t);
                end
            end
        end
        #if mod(sim,100) == 0 # for multiple sims keep output after every 100.
        #    println(sim)
        end
    return n,c,v
end


end



























































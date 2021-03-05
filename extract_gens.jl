module extract

export extract_gens, filter_lin

using Statistics

function extract_gens(dat_arr,col::Int)
    """ function that returns lineage data (2-d) split into the generations (3-d). """
    all_lins_data = [dat_arr[x][:,col] for x in 1:length(dat_arr)]; # extract the column of interest.
    div_flag_all_lins = [dat_arr[x][:,2] for x in 1:length(dat_arr)]; # the division flag data for 											each lineage.
    
    # now need to loop over all lineages and find the indices where division occurs.
    div_times = [findall(y-> y == 1.0, div_flag_all_lins[x]) for x in 1:length(div_flag_all_lins)];
    
    # init the 3-d generation sep data array.
    gen_sep_lins = [];
    
    # start loop over all lineages.    
    for x in 1:length(all_lins_data)
        gens = []; # init the 2-d generation measurements array.
        lin_data = all_lins_data[x]; # take the data for that lin.
        dt_lin = div_times[x]; # find the division times for that lin and then loop over.
        for y in 1:length(dt_lin)
            if y < length(dt_lin)
                push!(gens,lin_data[dt_lin[y]:dt_lin[y+1]-1]);
            else
                continue
            end
        end
        push!(gen_sep_lins,gens);
    end
    
    return gen_sep_lins
end


function filter_lin(lin::Array{Array{Float64,1},1}, n::Int64)
    """ function that applies an averaging filter over lin of width n """
    filter_lin = []; # init the filtered array.
    for (i,gen) in enumerate(lin)
        
        filter_gen = []; # the filter array for each gen.
        div = length(gen) ÷ n; # the int division of the gen length by the filter width
        rem = length(gen) % n; # the remainder

        for j in 1:div # loop over each full filter width
            j_prev = j-1;
            val = ones(n) .* mean(gen[j_prev*n+1:j*n])
            append!(filter_gen, val) # store to filter gen array
        end
        
        val_rem = ones(length(gen) - div*n) .* mean(gen[div*n+1:length(gen)]) # also include the remainder.
        
        append!(filter_gen, val_rem)
        push!(filter_lin, filter_gen) # push to the lineage filter array.
    end
    return convert(Array{Array{Float64,1},1},filter_lin), cat(filter_lin...,dims=1) # return filtered lin in sep and un-sep form.
end


end

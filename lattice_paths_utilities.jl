# Distributed for parallel stuff
# Delimited Files for csv export
using Distributed, DelimitedFiles
#addprocs(8) # uncomment to add slave processes

@everywhere begin
    using Random, SharedArrays, ProgressMeter

    function randomDyckPath(n) # uses cyclic lemma: build unbalanced path, cut at lowest point (with slight slope), rearrange and drop last item
        unb_walk = map(u->(u>n ? -1 : 1),randperm(2n+1))
        h = [1. * unb_walk[1]]
        for i in 1:length(unb_walk)
            push!(h,unb_walk[i]+h[end]+1/(2n+1))
        end
        return circshift(unb_walk,1-argmin(h))[1:end-1]
    end
    
    randomDyckBridge(n) = map(u->(u>n ? -1 : 1),randperm(2n))
    
    # pairs of well balanced steps/parenthesis
    function linksDyckPath(path)
        stack = []
        links = []
        h = 0
        for i in 1:length(path)
            p = path[i]                           # current step
            h += p                                # height of endpoint of current step
            sh = (h == 0 ? sign(h-p) : sign(h))   # if h positive, up steps open links and down steps close links.
                                                  # if h negative, viceversa
                                                  # if h zero, do as step before
            if path[i] == sh                      # if step opens link, psuh to stack
                push!(stack, i)
            else                                  # else pop and save link 
                j = pop!(stack)
                push!(links,[j,i])
            end
        end
        return links
    end
   
    # heights of up steps
    function heightsDyckPath(path) 
        heights = []
        h = 0
        for i in 1:length(path)
            step = path[i]
            h += step
            if step == 1
                push!(heights, h - step/2  )
            end
        end
        return heights
    end
   
    # entropy of path
    entropy(heights) = sum(heights .|> abs .|> u->u+0.5 |> log)
    
    # dyck cost of path (equispaced model)
    dyckCost(links,f=(u->u)) = foldl( (a,b)->a+f(b[2]-b[1]) , links, init=0 )   
    
    # transport field (equispaced model)
    transportfield(links) = map(b -> (b[2]-b[1]), sort(links,by=minimum) )

end

# prototype of simulation
#= function simulation(n,N,path) =#
#=     println("Start simulation \nNpts=",n,"\nNsim=",N) =#
#=     S = SharedArray{Float64}(N) =#
#=     @showprogress @distributed for i in 1:N =#
#=         S[i] = entropy(heightsDyckPath(randomDyckBridge(n))) =#
#=     end =#
#=     println("Done") =#
#=     output=replace(string( path,"entropy_brd_Npts",n,"_Nsim",N,".csv" )," "=>"") =#
#=     println("Saving at ",output) =#
#=     writedlm(output,S,',') =#
#=     println("Done") =#
#= end =#
#= for n in [1000000] =#
#=     @time simulation(n,100000,"./data/") =#
#=     println(repeat("#",40)) =#
#= end =#

# ASEN 6519 Final Project
# Collin Hudson 12/07/2024
using Plots
using Printf
using Distributions
using Zygote
using LazySets
using JuMP, SCIP
using TickTock
ENV["TICKTOCK_MESSAGES"] = false
if isdefined(@__MODULE__, :workspace)
    workspace isa Module || error("workspace is present and it is not a Module")
else
    include("workspace.jl")
end
using .workspace
# Workspace parameters (xlim, ylim, number of obstacles)
    # seed 2300 = squeeze
    # 353 broken somehow
    e = env([-5,5],[-5,5],6,62)
    # Nverts = sum(x->length(x.vertices),e.obs)
# 

# # Solve for optimal M values for obstacle avoidance
#     ks = zeros(Nverts,1)
#     hs = zeros(Nverts,2*Nverts)
#     row = 1
#     for ob in e.obs
#         for h in LazySets.constraints_list(ob)
#             hs[row,:] = hcat(zeros(1,2*(row-1)),h.a',zeros(1,2*(Nverts - row)))
#             ks[row] = h.b
#             global row += 1
#         end
#     end
#     model = Model(HiGHS.Optimizer)
#     set_silent(model)
#     @variable(model, x[i=1:2*Nverts])
#     for j in 1:2:2*Nverts
#         @constraint(model,e.xlim[1] <= x[j] <= e.xlim[2])
#         @constraint(model,e.ylim[1] <= x[j+1] <= e.ylim[2])
#     end
#     @objective(model, Max, sum(ks - hs*x))
#     optimize!(model)
#     # Mobs = maximum(ks - hs*value.(x))*ones(length(value.(x)),1)
#     Mobs = ks - hs*value.(x)
# # 

# Problem parameters
    Nga = 3
    maxUAS = 5
    Ntot = Nga+maxUAS
    rCon = 3
    rMin = 0
    rUAS = 2
    # # gaPosX = [[1;0;3],[0;1.5;2],[-0.5;1.5;0]]
    # # gaPosY = [[2;-2;-4],[2;-1.5;-4],[3;0;-4]]
    gaPosX = [[1;0;5],[0;-1;4],[0;-2;3]]
    gaPosY = [[2;-1.5;-1],[2.5;-1.5;-1],[3;-1.8;-2]]
    samples = 0.5
    # Mcon = 1.1*((10)^2 + (10)^2)
    # # Mcon = 1.1*sqrt((e.xlim[2] - e.xlim[1])^2 + (e.ylim[2] - e.ylim[1])^2)
    Mcon = 1.1*((e.xlim[2] - e.xlim[1])^2 + (e.ylim[2] - e.ylim[1])^2)
#

function findStep(e, xLast, yLast, xcLast, k)
    model = Model(SCIP.Optimizer)
    set_silent(model)
    set_time_limit_sec(model,45)
    # Variables and workspace boundary constraints
        @variable(model, e.xlim[1] <= x[i=1:Ntot] <= e.xlim[2]) #Agent positions x
        @variable(model, e.ylim[1] <= y[i=1:Ntot] <= e.ylim[2]) #Agent positions y
        @variable(model, c[i=1:Ntot,j=1:Ntot], Bin) #Network edge choice (which nodes are connected)
        @variable(model, o[i=1:e.Nverts,j=1:Ntot],Bin) #Slack variables for obstacle avoidance
        @variable(model, l[r=1:e.Nverts,i=1:Ntot,j=1:Ntot,k=1:length(samples)],Bin) #Slack variables for LOS
        @variable(model, xc[i=1:Ntot],Bin) #Network cluster (active connection set)
        @variable(model, f[i=1:Ntot,j=1:Ntot] >= 0,Int) #Network flow
        @objective(model, Min, sum(xc) + sum((xc - xcLast).^2) + 0.1*sum(c[:,(Nga+1):Ntot]*ones(maxUAS,1))) #Minimize number of agents in network cluster and number of connections to UAS
    #

    # Ground agent constraints
        for j in 1:Nga
            fix(x[j], xLast[j];force=true) #Ground agents must be at given x
            fix(y[j], yLast[j];force=true) #Ground agents must be at given y
            fix(xc[j], 1;force=true) #Ground agents must be in network cluster
        end
    # 

    # Timestep constraints
        if k != 1
            for i in (Nga+1):Ntot
                @constraint(model, (x[i] - xLast[i])^2 + (y[i] - yLast[i])^2 <= (rUAS^2 + Mcon*(1-xcLast[i])))
                # fix(xc[i], xcLast[i];force=true)
            end
            # Add LOS constraint by adding columns? to o for each last position uas
        end
    # 

    # Network constraints
        @constraint(model, sum(c*ones(Ntot,1)) >= (sum(xc)-1)) #Network must have at least (number of cluster agents - 1) edges
        for i in 1:Ntot
            @constraint(model, sum(c[i,:]) <= Ntot*xc[i]) #Only agents in cluster can be connected
            @constraint(model, sum(c[:,i]) <= Ntot*xc[i]) #Only agents in cluster can be connected
            @constraint(model, sum(f[i,:]) <= Ntot*xc[i]) #Only agents in cluster can have flow
            # @constraint(model, sum(f[:,i]) <= Ntot*xc[i]) #Only agents in cluster can have flow
            if i != 1
                @constraint(model,(sum(f[i,:]) - sum(f[:,i])) == -1*xc[i]) #Non-source cluster agents must consume one unit of flow
            end
            for j in 1:Ntot
                if i < j
                    @constraint(model, (f[i,j] + f[j,i]) <= ((sum(xc) - 1)*c[i,j]))
                elseif i == j
                    fix(f[i,j], 0;force=true) #Flow must be between two different agents
                end
            end
        end
    # 

    # Connection radius constraint
        for i in 1:Ntot
            for j in 1:Ntot
                if i < j
                    @constraint(model, (x[i] - x[j])^2 + (y[i] - y[j])^2 <= (rCon^2 + Mcon*(1-c[i,j]) + Mcon*(1-xc[j])))
                    @constraint(model, (x[i] - x[j])^2 + (y[i] - y[j])^2 >= xc[j]*rMin^2)
                else
                    fix(c[i,j], 0;force=true)
                end
            end
        end
    # 

    # Obstacle avoidance and LOS constraints
        row = 1
        obV = 1
        for ob in e.obs
            for con in constraints_list(ob)
                for j in 1:Ntot
                    @constraint(model,-con.a'*[x[j];y[j]] <= (-con.b + e.Mobs[row]*o[row,j]))
                end
                for i in 1:Ntot
                    for j in 1:Ntot
                        if i < j
                            for k in eachindex(samples)
                                si = (1-samples[k])*[x[i];y[i]] + samples[k]*[x[j];y[j]]
                                @constraint(model, -con.a'*si <= (-con.b + e.Mobs[row]*o[row,i] + e.Mobs[row]*l[row,i,j,k]))
                                @constraint(model, -con.a'*si <= (-con.b + e.Mobs[row]*o[row,j] + e.Mobs[row]*l[row,i,j,k]))
                            end
                            @constraint(model, c[i,j] <= sum(1 .- l[row,i,j,1:length(samples)])) #At least one sample point must be valid if connected
                        end
                    end
                end
                row += 1
            end
            for j in 1:Ntot
                @constraint(model,sum(o[obV:(obV+length(ob.vertices)-1),j]) <= (length(ob.vertices)-1)) #At least one constraint must be active per obstacle
            end
            obV += length(ob.vertices)
        end
    # 

    optimize!(model)
    if !is_solved_and_feasible(model)
        p = plot(size = (400,400))
        for j in eachindex(e.obs)
            plot!(p,e.obs[j])
        end
        xlims!(p,e.xlim[1]-1,e.xlim[2]+1)
        ylims!(p,e.ylim[1]-1,e.ylim[2]+1)
        for j in 1:Ntot
            if j < Nga+1
                scatter!([xLast[j]],[yLast[j]],mc=:green,label=nothing)
                plot!(xLast[j] .+ rCon*cos.(range(0,2*π,500)),yLast[j] .+ rCon*sin.(range(0,2*π,500)),linecolor=:green,label=nothing)
            else
                scatter!([xLast[j]],[yLast[j]],mc=:black,label=nothing)
                plot!(xLast[j] .+ rUAS*cos.(range(0,2*π,500)),yLast[j] .+ rUAS*sin.(range(0,2*π,500)),linecolor=:gray,label=nothing)
            end
        end
        title!(p,"Infeasible workspace")
        display(p);
    end
    if is_solved_and_feasible(model)
        return [value.(x),value.(y),ceil.(value.(xc)),ceil.(value.(c))]
    else
        return [zeros(1,Ntot),zeros(1,Ntot),zeros(1,Ntot),zeros(Ntot,Ntot)]
    end
end

# solve loop
    # solve for first timestep
    xPath = []
    yPath = []
    xcPath = []
    cPath = []
    for k in eachindex(gaPosX)
        if k == 1
            (xk,yk,xck,ck) = findStep(e,[gaPosX[1];zeros(maxUAS,1)],[gaPosY[1];zeros(maxUAS,1)],zeros(Ntot,1),1)
        else
            (xk,yk,xck,ck) = findStep(e,[gaPosX[k];xPath[k-1][(Nga+1):Ntot]],[gaPosY[k];yPath[k-1][(Nga+1):Ntot]],xcPath[k-1],k)
        end
        push!(xPath,xk)
        push!(yPath,yk)
        push!(xcPath,xck)
        push!(cPath,ck)
    end

#  Monte Carlo
    # N = 50
    # nobsV = [2,3,4]
    # xcMonte = zeros(N,length(nobsV))
    # timeMonte = zeros(N,length(nobsV))
    # for i in eachindex(nobsV)
    #     idx = 1
    #     tries = 0
    #     global Nga = nobsV[i]
    #     global Ntot = Nga+maxUAS
    #     while (idx < (N+1)) && tries < 100
    #         # Generate random environment
    #         global e = env([-5,5],[-5,5],5)
    #         badGA = true
    #         triesEnv = 0
    #         while badGA && triesEnv < 50
    #             global gaPosX = 2*e.xlim[2] .*rand(nobsV[i]) .+ e.xlim[1]
    #             global gaPosY = 2*e.ylim[2] .*rand(nobsV[i]) .+ e.ylim[1]
    #             badGA = false
    #             for o in e.obs
    #                 for j in 1:Nga
    #                     badGA = badGA||([gaPosX[j];gaPosY[j]] ∈ o)
    #                 end
    #             end
    #             triesEnv += 1
    #         end
    #         tick()
    #         if !badGA
    #             (xk,yk,xck,ck) = findStep(e,[gaPosX;zeros(maxUAS,1)],[gaPosY;zeros(maxUAS,1)],zeros(Ntot,1),1)
    #             if sum(xck) >= Nga
    #                 if sum(xck) > Ntot
    #                     print(xck)
    #                 end
    #                 xcMonte[idx,i] = sum(xck) - Nga
    #                 timeMonte[idx,i] = tok()
    #             else
    #                 xcMonte[idx,i] = -99
    #                 timeMonte[idx,i] = -tok()
    #             end
    #             idx += 1
    #         else
    #             tok()
    #         end
    #         print("Tries: "*string(tries)*", monte: "*string(idx-1)*" Xga: "*string(gaPosX')*" Yga: "*string(gaPosY')*"\n")
    #         tries += 1
    #     end
    # end

    
    
# 


# # Plotting
    for k in eachindex(xPath)
        p = plot(size = (400,400))
        for j in eachindex(e.obs)
            plot!(p,e.obs[j])
        end
        xSol = xPath[k]
        ySol = yPath[k]
        cSol = cPath[k]
        xcSol = xcPath[k]
        for j in 1:Ntot
            for i in 1:Ntot
                if(cSol[j,i] != 0) && xcSol[j] != 0
                    plot!([xSol[j],xSol[i]], [ySol[j],ySol[i]],linecolor=:gray,label=nothing)
                end
            end
            if j <= Nga
                scatter!([xSol[j]],[ySol[j]],mc=:green,label=nothing)
            else
                if xcSol[j] != 0
                    scatter!([xSol[j]],[ySol[j]],mc=:yellow,label=nothing)
                    # plot!(xSol[j] .+ rCon*cos.(range(0,2*π,500)),ySol[j] .+ rCon*sin.(range(0,2*π,500)),linecolor=:gray,label=nothing)
                else
                    scatter!([xSol[j]],[ySol[j]],mc=:black,label=nothing)
                end
            end
        end
        xlims!(p,e.xlim[1]-1,e.xlim[2]+1)
        ylims!(p,e.ylim[1]-1,e.ylim[2]+1)
        for j in 1:Nga
            plot!(gaPosX[k][j] .+ rCon*cos.(range(0,2*π,500)),gaPosY[k][j] .+ rCon*sin.(range(0,2*π,500)),linecolor=:green,label=nothing)
        end
        title!(p,"Timestep: "*string(k))
        display(p);
    end
# # 

# Author: Noé Blassel
# Date: 2025/06

using Distributed

println("Hello from main process at $(gethostname()) !")
machines = readlines(ENV["OAR_NODEFILE"])
hostname = gethostname()
ssh_agent = "oarsh"

println("$(length(machines)) available cores across $(length(Set(machines))) nodes")

flush(stdout)
addprocs(machines,ssh=ssh_agent)

@everywhere basedir = pwd()

if !isdir("trajectories")
    mkdir(joinpath(basedir,"trajectories"))
end

@everywhere function umbrella_run(ϕ,ψ,nsteps)
    files = ["aladim_amber.xyz","aladim_amber.key","amber99.prm","amber99sb.prm","aladim_amber.tcl"]
    junk_files = ["aladim_amber.colvars.state","aladim_amber.colvars.state.old","aladim_amber.colvars.traj.BAK"]

    workingdir = joinpath(basedir,"umbrella$(ϕ)_$(ψ)")

    !isdir(workingdir) && mkdir(workingdir) # make local directory for output
    cd(workingdir)

    for f in files
        cp(joinpath(basedir,f),joinpath(workingdir,f),force=true)
    end

    dtsteer = dtrun = 1.0
    steered = success = false
    nsteer = 750
    nrun = nsteps

    nfailsteer = nfailrun = 0
    max_try = 5
    
    while !steered || !success
        ## try steering run
        try
            f = joinpath(basedir,"templates","steer.colvars")
            s = read(f,String)
            w = open("aladim_amber.colvars","w")
            print(w,replace(s,"XXX YYY"=>"$ϕ $ψ","NNN"=>"$nsteer")) # set harmonic center
            close(w)
            run(`/home/nblassel/tinker-hp/v1.2/bin/dynamic aladim_amber $nsteer $dtsteer 1 2 300`) # run steered md
            steered = true
        catch # if steering run fails
            dtsteer /= 2
            nsteer *= 2

            for f in junk_files
                rm(f,force=true)
            end

            rm("aladim_amber.dyn",force=true)
            rm("aladim_amber.dcd",force=true)

            nfailsteer += 1
        end

        for f in junk_files
            rm(f,force=true)
        end

    # if steering successful, move on to sampling run
        if steered
            try
                f = joinpath(basedir,"templates","harm.colvars")
                s = read(f,String)
                w = open("aladim_amber.colvars","w")
                print(w,replace(s,"XXX YYY"=>"$ϕ $ψ"))
                close(w)
                run(`/home/nblassel/tinker-hp/v1.2/bin/dynamic aladim_amber $nrun $dtrun 1 2 300`) # sampling run
                success = true
            catch # if sampling run fails
                dtrun /= 2
                nrun *= 2

                for f in junk_files
                    rm(f,force=true)
                end

                rm("aladim_amber.dyn",force=true)
                rm("aladim_amber.dcd",force=true)

                steered = false # fall back to steered md run
                nfailrun += 1
            end

            if ((nfailsteer > max_try) || (nfailrun > max_try))
                println("Error: steered MD failed $nfailsteer times and umbrella sampling failed $nfailrun times. Moving on.")
                break
            end
        end
    end

    try
        cp(joinpath(workingdir,"aladim_amber.colvars.traj"),joinpath(basedir,"trajectories","traj_$(ϕ)_$(ψ).out"),force=true)

        cd(basedir)
        rm(workingdir,force=true,recursive=true) # cleanup
    catch
    end
end

ϕrange = -235:5:125
ψrange = -120:5:240

centers = [(ϕ,ψ) for ϕ=ϕrange,ψ=ψrange]
println("$(length(centers)) umbrella sampling runs to complete over $(nworkers()) workers.")
nsteps = 50000 # EDIT HERE

pmap(t->umbrella_run(t...,nsteps),centers;on_error = e->println(e)) # distribute umbrella sampling runs
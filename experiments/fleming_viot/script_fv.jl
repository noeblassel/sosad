#!/home/nblassel/.juliaup/bin/julia

# replace with path

using DelimitedFiles

# Define phases
@enum Phase Decorr Sample

# Parse arguments
if length(ARGS) < 7
    println(stderr, "Usage: script.jl <isaddle> <nreps> <nsteps_decorr> <nsteps_times> <isopt> <Nruns> <gamma>")
    exit(1)
end
isaddle, nreps, nsteps_decorr, nsteps_times, isopt, Nruns, gamma = parse.(Int, ARGS)

# Constants and directories
dt = 2.0
saddles = Dict(1=>(-3,89), 2=>(44,-92), 3=>(61,122), 4=>(124,-19),5=>(53,28)) # saddle 5 = free-energy minimum
ϕm,ψm = saddles[5]
ϕ, ψ = saddles[isaddle]
suffix = isopt == 0 ? "_basin" : "_opt"

basedir = pwd()
templates = joinpath(basedir, "templates")
initdir = mkpath(joinpath(basedir, "init$(suffix)_$(gamma)","saddle$(ϕ)_$(ψ)"))

tclscript = "aladim_amber_states.tcl"
cv_file_run = "fv"*suffix*".colvars"


# reperiodize to nearest periodic image wrt free-energy minimum
reperiodize_phi(ϕ) = begin
    if abs(ϕ-ϕm-360) < abs(ϕ-ϕm)
        return ϕ-360
    elseif abs(ϕ-ϕm+360) < abs(ϕ-ϕm)
        return ϕ+360
    else
        return ϕ
    end
end

reperiodize_psi(ψ) = begin
    if abs(ψ-ψm-360) < abs(ψ-ψm)
        return ψ-360
    elseif abs(ψ-ψm+360) < abs(ψ-ψm)
        return ψ+360
    else
        return ψ
    end
end

# Steered MD
function steered_md(ϕ, ψ; steps=2000, dt=0.5)
    try
        cfg = read(joinpath(templates, "steer.colvars"), String)
        cfg = replace(cfg, "XXX YYY"=>"$ϕ $ψ", "NNN"=>"$steps")
        write("aladim_amber.colvars", cfg)
        run(`/home/nblassel/tinkerhp-genparrep/CPU/bin/dynamic aladim_amber $steps $dt 1 2 300`) # replace with tinkerhp-genparrep
        return true
    catch e
        println(stderr, "[Error] Steered MD failed: $e")
        flush(stderr)
        return false
    end
end

# Harmonic relaxation MD
function harmonic_md(ϕ, ψ; steps=5000,checksteps=2000, dt=1.0)
    try
        check_flag = false
        firsttry_flag = true

        while !check_flag
            clean_dir()
            cfg = read(joinpath(templates, "harm.colvars"), String)
            cfg = replace(cfg, "XXX YYY"=>"$ϕ $ψ")
            write("aladim_amber.colvars", cfg)

            if firsttry_flag
                run(`/home/nblassel/tinkerhp-genparrep/CPU/bin/dynamic aladim_amber $steps $dt 1 2 300`)
                firsttry_flag = false
            else
                run(`/home/nblassel/tinkerhp-genparrep/CPU/bin/dynamic aladim_amber $checksteps $dt 1 2 300`)
            end

            check_flag = check() # check system has entering velocity
        end
        return true
    catch e
        println(stderr, "[Error] Harmonic MD failed: $e")
        flush(stderr)
        return false
    end
end

# Check system has entering velocities and is in free-energy basin
function check()
    (isaddle == 5) && return true
    try

        cp(joinpath(templates,"check.colvars"),"aladim_amber.colvars";force=true)
        cp(joinpath(basedir,tclscript),tclscript;force=true)

        for f in ["aladim_amber.colvars.traj","aladim_amber.colvars.traj.BAK","aladim_amber.colvars.state"]
            rm(f;force=true)
        end

        run(`/home/nblassel/tinkerhp-genparrep/CPU/bin/dynamic aladim_amber 2 1 0.001 2 300`) # run check -- two steps of equilibrium Langevin

        M = readdlm("aladim_amber.colvars.traj";comments=true)
        ϕ1,ϕ2 = reperiodize_phi.(M[1:2,2])
        ψ1,ψ2 = reperiodize_psi.(M[1:2,3])

        statecheck = (M[2,4]>0)
        velcheck = (sum(abs2,[ϕ1-ϕm,ψ1-ψm]) > sum(abs2,[ϕ2-ϕm,ψ2-ψm]))

        !velcheck && println(stderr, "[Warning] in prep check: velocity is exiting. Proceeding to next iteration.")
        !statecheck && println(stderr, "[Warning] in prep check: initial point was not in basin ($(M[2,:])). Proceeding to next iteration.")

        return   velcheck && statecheck

    catch e
        println(stderr, "[Error] in prep check: $e")
        flush(stderr)
        return false
    end
end

# Prepare initial condition for a trial
tmp_init(k) = joinpath(initdir, "init_$(isaddle)_trial$(k).dyn")

function prepare_initial_trial(k)
    prep_dir = mkpath(joinpath(basedir, "prep$(ϕ)_$(ψ)$(suffix)_trial$(k)_$(gamma)"))
    success = false
    cd(prep_dir) do
        for f in ("aladim_amber.xyz", "amber99.prm", "amber99sb.prm","aladim_amber.key")
            cp(joinpath(basedir, f), f; force=true)
        end
        if steered_md(ϕ, ψ) && harmonic_md(ϕ, ψ)
            cp("aladim_amber.dyn", tmp_init(k); force=true)
            success = true
        else
            println(stderr, "[Error] Initial prep failed at trial $k.")
        end
    end
    rm(prep_dir; force=true, recursive=true)
    return success
end

# Read trajectories
function read_traj(j)
    idx = lpad(string(j), 3, '0')
    fname = "aladim_amber_reps" * idx * ".colvars.traj"
    data = readdlm(fname; comments=true)
    return data[:,2:end]
end

function clean_dir()
    try
        run(`sh -c "rm *.dcd"`)
        run(`sh -c "rm *.state*"`)
    catch e
    end
end

# Run one phase
function run_phase(phase::Phase, steps, k)
    try
        dir = joinpath(basedir, "trajectories$(suffix)_$(gamma)", "saddle"*string(ϕ)*"_"*string(ψ))
        mkpath(dir)
        j = length(readdir(dir)) ÷ 4 # to avoid overwritting data from previous runs

        if j==Nruns
            println("Master script: target number of runs achieved.")
            exit(0)
        end

        run(`mpirun --mca orte_rsh_agent "oarsh" -machinefile $(get(ENV,"OAR_NODEFILE","")) -n $nreps -x LD_LIBRARY_PATH -x PATH /home/nblassel/tinkerhp-genparrep/CPU/bin/dynamic aladim_amber $steps $dt 1 2 300`)
        data = [read_traj(j) for j in 0:nreps-1]


        out = hcat(data...)

        if phase == Decorr
            writedlm(joinpath(dir, "decorr_fv_$(j+1).out"), out)
            cp("fv_exit_times.out",joinpath(dir,"times_decorr_$(j+1).out");force=true)

            run(`sh -c "rm *.colvars.state*"`)
        else
            writedlm(joinpath(dir, "eq_fv_$(j+1).out"), out)
            cp("fv_exit_times.out", joinpath(dir, "times_eq_$(j+1).out"); force=true)

            clean_dir()
        end

        return true
    catch e
        start_cond = tmp_init(k)
        println(stderr, "[Error] Phase=$(phase) trial=$(k) failed: $(e). Deleting $start_cond")
        flush(stderr)
        return false
    end
end

# Run a full FV trial
function run_fv_trial(k)
    # sample new initial condition
    if !prepare_initial_trial(k)
        return false
    end
    trial_dir = mkpath(joinpath(basedir,"fv_runs$(suffix)_$(gamma)", "saddle$(ϕ)_$(ψ)", "trial_$(k)"))

    ok1 = ok2 = false

    cd(trial_dir) do
        for f in ("aladim_amber.xyz", "amber99.prm", "amber99sb.prm",tclscript)
            cp(joinpath(basedir, f), f; force=true)
        end
        cp(tmp_init(k), "aladim_amber.dyn"; force=true)
        key_template = read(joinpath(templates, "fv.key"), String)
        write("aladim_amber.key", replace(key_template, "NNN"=>string(nreps),"FFF"=>string(float(gamma))))
        cp(joinpath(templates,cv_file_run),"aladim_amber.colvars"; force=true)
        ok1 = run_phase(Decorr, nsteps_decorr, k)
        ok2 = ok1 && run_phase(Sample, nsteps_times, k)
    end

    rm(trial_dir;force=true,recursive=true)
    return ok1 && ok2
end

function main()
    # Main workflow
    println("=== Starting FV sampling with dynamic inits ===")
    success = 0
    attempt = 1
    logfile = open("fv_run$(suffix)_$(ϕ)_$(ψ)_$(gamma).log","w")
    while success < Nruns
        println("Attempt $(attempt) Successes=$(success)/$(Nruns)")
        if run_fv_trial(attempt)
            success += 1
        end
        attempt += 1
        println(logfile,"Completed $(success) successful runs in $(attempt-1) attempts.")
        flush(logfile)
    end
end

main()
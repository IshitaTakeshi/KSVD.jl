# The code is from
# http://thirld.com/blog/2015/05/30/julia-profiling-cheat-sheet/
function benchmark(func, args...; filename = "trace.def", delay = 0.01)
    # Any setup code goes here.

    # Run once, to force compilation.
    println("======================= First run:")
    @time func(args...)

    # Run a second time, with profiling.
    println("======================= Second run:")

    Profile.init(delay = delay)
    Profile.clear()
    Profile.clear_malloc_data()

    @profile @time func(args...)

    # Write profile results to profile.bin.
    f = open(filename, "w")
    Profile.print(f)
    close(f)
end

using Random, Base.Threads

# Parallel random number generation
n = 10000
random_nums = Vector{Float64}(undef, n)

@threads for i in 1:n
    random_nums[i] = rand()
end

println(random_nums)

# using Random, Base.Threads

# n = 5
# random_nums = Vector{Float64}(undef, n)

# # Seed each thread's RNG
# @threads for i in 1:n
#     rng = MersenneTwister(i)  # Seeded RNG per thread
#     random_nums[i] = rand(rng)
# end

# println(random_nums)


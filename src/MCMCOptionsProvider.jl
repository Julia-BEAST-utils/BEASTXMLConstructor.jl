module MCMCOptionsProvider

export MCMCOptions

mutable struct MCMCOptions
    chain_length::Int
    file_log_every::Int
    screen_log_every::Int
    likelihood_check_count::Int
end

function MCMCOptions(;
                     chain_length::Int = 10_000,
                     file_log_every::Int = max(1, div(chain_length, 1_000)),
                     screen_log_every::Int = max(1, div(chain_length, 100)),
                     likelihood_check_count::Int = -1)

    if likelihood_check_count == -1
        check_count = div(chain_length, 10)
        check_count = max(check_count, 100)
        check_count = min(check_count, 1000)
    else
        check_count = likelihood_check_count
    end

    return MCMCOptions(chain_length, file_log_every, screen_log_every,
                       likelihood_check_count)
end

end
"""
    runmean = running_mean(v, w)

return the running mean (`runmean`) of `v` using kernel width `w`
"""
function running_mean(v, w::Number) 
    runmean = [i < w ? Statistics.mean(v[begin:i]) : Statistics.mean(v[i-w+1:i]) for i in 1:length(v)]
    return runmean
end
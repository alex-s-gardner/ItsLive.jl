"""
List all .nc files under s3://its-live-data-test/velocity_image_pair/nisar/v02/
Uses AWSS3.jl + AWS.jl with anonymous access (public bucket, no credentials needed).
"""

using AWS, AWSS3

const NISAR_BUCKET = "its-live-data-test"
const NISAR_PREFIX = "velocity_image_pair/nisar/v02/"

function list_nisar_nc_files(
    bucket::String = NISAR_BUCKET,
    prefix::String = NISAR_PREFIX;
    region::String = "us-west-2",
    verbose::Bool = false,
)
    cfg = AWSConfig(; creds=nothing, region=region)
    verbose && println("Listing s3://$bucket/$prefix ...")

    nc_files = String[]
    for obj in s3_list_objects(cfg, bucket, prefix; delimiter="")
        key = obj["Key"]
        endswith(key, ".nc") && push!(nc_files, "s3://$bucket/$key")
    end

    return nc_files
end

"""
    s3_to_https(s3_url; region="us-west-2") -> String

Convert an S3 URL (`s3://bucket/key`) to an HTTPS URL for the AWS S3 REST endpoint.

# Example
```julia
s3_to_https("s3://its-live-data-test/velocity_image_pair/nisar/v02/N20E080/foo.nc")
# → "https://its-live-data-test.s3.us-west-2.amazonaws.com/velocity_image_pair/nisar/v02/N20E080/foo.nc"
```
"""
function s3_to_https(s3_url::String; region::String = "us-west-2")
    startswith(s3_url, "s3://") || error("Not an s3:// URL: $s3_url")
    rest   = s3_url[6:end]                          # strip "s3://"
    slash  = findfirst('/', rest)
    isnothing(slash) && error("No key found in S3 URL: $s3_url")
    bucket = rest[1:slash-1]
    key    = rest[slash+1:end]
    return "https://$bucket.s3.$region.amazonaws.com/$key"
end

# Run when executed directly
files = list_nisar_nc_files(; verbose=true)
println("\nFound $(length(files)) .nc files:\n")
foreach(println, files)

path = s3_to_https.(files)


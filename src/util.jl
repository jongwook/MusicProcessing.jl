# utility functions

tofloat(array::AbstractArray{T}) where T = map(Float32, array)
tofloat(array::AbstractArray{Float32}) = array

fft_frequencies(samplerate::Real, nfft::Int) = collect(range(0f0, samplerate / 2f0, length=(nfft >> 1) + 1))


"""returns the number of frames when the signal is partitioned into overlapping frames"""
nframes(length::Int, framesize::Int, hopsize::Int) = div(length - framesize, hopsize) + 1

"""Provides a view of a signal partitioned into overlapping frames, similar to MATLAB's buffer"""
struct FrameView{T<:AbstractVector} <: AbstractVector{Vector}
    original::T
    framesize::Int
    hopsize::Int
    nframes::Int

    function FrameView(original, framesize, hopsize)
        # n = noverlap is a problem - the algorithm will not terminate.
        #@boundscheck((0 < hopsize <= framesize), error("hopsize must be between zero and framesize"))
        new{T}(original, framesize, hopsize, nframes(length(original), framesize, hopsize))
    end
end

function Base.getindex(x::FrameView, i::Int)
    #@boundscheck((i >= 1 && i <= x.k), Base.throw_boundserror(x, i))
    slice(x.original, (i-1) * x.hopsize + (1:x.framesize))
end

Base.size(x::FrameView) = (x.nframes,)

function Base.iterate(x::FrameView, i::Int)
   if i > x.nframes
       nothing
   else
       (x[i], i+1)
   end
end



"""returns the DCT filters"""
function dct(nfilters::Int, ninput::Int)
    basis = Array{Float32}(undef, nfilters, ninput)
    samples = (1f0:2f0:2ninput) * π / 2ninput
    for i = 1:nfilters
        for j in 1:length(samples)
            basis[i, j] = cos(i * samples[j])
        end
    end

    basis *= sqrt(2f0/ninput)
    basis
end

"""returns the number of times that the value changes its sign"""
function zero_crossings(array::AbstractVector{T}, size::Int = length(array), offset::Int = 0) where {T<:Real}
    result = 0
    previous = 0
    for i in offset .+ (1:size)
        number = array[i]
        sgn = number == 0 ? 0 : number > 0 ? 1 : -1
        if sgn != previous
            result += 1
        end
        previous = sgn
    end
    result
end

# units from SampledSignals

const frames = Hz*s
const FrameQuant = DimensionlessQuantity

"""
    inframes([Type,]quantity[, rate])

Translate the given quantity to a (unitless) number of time or frequency frames,
given a particular samplerate. Note that this isn't quantized to integer numbers
of frames. If given a `Type`, the result will first be coerced to the given type.

If the given quantity is Unitful, we use the given units. If it is not we assume
it is already a value in frames.

# Example

julia> inframes(0.5s, 44100Hz)
22050.0

julia> inframes(1000Hz, 2048/44100Hz)
46.439909297052154

"""
inframes(::Type{T}, frame::FrameQuant, rate=nothing) where T <: Integer =
    round(T, ustrip(uconvert(frames, frame)))
inframes(frame::FrameQuant, rate=nothing) = ustrip(uconvert(frames, frame))
inframes(::Type{T}, time::Unitful.Time, rate) where T <: Integer =
    round(T, inseconds(time)*inHz(rate))
inframes(time::Unitful.Time, rate) = inseconds(time)*inHz(rate)
inframes(::Type{T}, freq::Unitful.Frequency, rate) where T <: Integer =
    round(T, inHz(freq)*inseconds(rate))
inframes(freq::Unitful.Frequency, rate) = inHz(freq)*inseconds(rate)
inframes(::Type, frame::Unitful.AbstractQuantity) = error("Unknown sample rate")
inframes(frame::Unitful.AbstractQuantity) = error("Unknown sample rate")
inframes(::Type{T}, frame::Real) where T = T(frame)
inframes(frame::Real) = frame

"""
    inHz(quantity[, rate])

Translate a particular quantity (usually a frequency) to a (unitless) value in
Hz.

If the given quantity is Unitful, we use the given units. If it is not we assume
it is already a value in Hz.

For some units (e.g. frames) you will need to specify a sample rate:

# Example

julia> inHz(1.0kHz)
1000.0

"""
inHz(x::Unitful.Frequency, rate=nothing) = ustrip(uconvert(Hz, x))
inHz(x::FrameQuant) = error("Unknown sample rate")
# assume we have a spectrum buffer with a sample rate in seconds
inHz(x::FrameQuant, rate) = inHz(inframes(x) / rate)
inHz(x::Real, rate) = x
inHz(x::Real) = x

"""
   inseconds(quantity[, rate])

Translate a particular quantity (usually a time) to a (unitless) value in
seconds.

If the given quantity is Unitful, we use the given units. If it is not we assume
it is already a value in seconds.

For some units (e.g. frames) you will need to specify a sample rate:

# Examples
julia> inseconds(50.0ms)
0.05

julia> inseconds(441frames, 44100Hz)
0.01

"""
inseconds(x::Unitful.Time, rate=nothing) = ustrip(uconvert(s,x))
inseconds(x::FrameQuant) = error("Unknown sample rate")
# assume we have a time buffer with sample rate in hz
inseconds(x::FrameQuant, rate) = inseconds(inframes(x) / rate)
inseconds(x::Real, rate) = x
inseconds(x::Real) = x

# inseconds(x, rate::Quantity) = inseconds(x,inHz(rate))
# inseconds(x::FrameQuant, rate::Real) = (ustrip(x) / rate)
# inseconds(x::Quantity, rate::Real) = ustrip(uconvert(s,x))

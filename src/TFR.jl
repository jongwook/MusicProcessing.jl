
import Base.eps
import DSP.spectrogram, DSP.stft

""""""
function spectrogram(audio::SampleBuf{T, 1},
                        windowsize::Int = 1024,
                        hopsize::Int = windowsize >> 2;
                        window = hanning, kwargs...) where T
    noverlap = windowsize - hopsize
    DSP.spectrogram(audio.data, windowsize, noverlap; fs = audio.samplerate, window = window, kwargs...)
end

""""""
function spectrogram(audio::SampleBuf{T, 2},
                        windowsize::Int = 1024,
                        hopsize::Int = windowsize >> 2;
                        window = hanning, kwargs...) where T
    noverlap = windowsize - hopsize
    vec(mapslices(audio.data, dims=1) do data
        DSP.spectrogram(data, windowsize, noverlap; fs = audio.samplerate, window = window, kwargs...)
    end)
end

""""""
function spectrogram(audio::SampleBuf{T, N},
                           windowsize::Seconds,
                           hopsize::Seconds = windowsize / 4;
                           kwargs...) where {T, N}
    w = round(Int, ustrip(windowsize) * audio.samplerate)
    h = round(Int, ustrip(hopsize) * audio.samplerate)
    spectrogram(audio, w, h; kwargs...)
end


""""""
function stft(audio::SampleBuf{T, 1},
                 windowsize::Int = 1024,
                 hopsize::Int = windowsize >> 2;
                 window = hanning, kwargs...) where T
    noverlap = windowsize - hopsize
    DSP.stft(audio.data, windowsize, noverlap; window = window, kwargs...)
end

""""""
function stft(audio::SampleBuf{T, 2},
                 windowsize::Int = 1024,
                 hopsize::Int = windowsize >> 2; kwargs...) where T
    nchannels = size(audio.data, 2)
    noverlap = windowsize - hopsize

    stft = Array{Matrix{Complex{Float32}}}(undef, nchannels) # type that DSP.stft outputs
    for i in 1:nchannels
        stft[i] = DSP.stft(audio.data[:, i], windowsize, noverlap; kwargs...)
    end
    cat(stft..., dims=3)
end

""""""
function stft(audio::SampleBuf{T, N},
                    windowsize::Seconds,
                    hopsize::Seconds = windowsize / 4; kwargs...) where {T, N}
    w = round(Int, ustrip(windowsize) * audio.samplerate)
    h = round(Int, ustrip(hopsize) * audio.samplerate)
    stft(audio, w, h; kwargs...)
end


# functions to make window vector from a window function, vector or nothing.
@inline w(window::Function, size) = window(size)
@inline function w(window::AbstractVector, size)
    if size != length(window)
        error("the length of window does not match windowsize")
    end
    window
end
@inline w(window::Nothing, size) = ones(size)


"""
    istft(stft, samplerate, windowsize, hopsize; window)

Reconstruct a mono audio from its STFT.

# Arguments
* `stft::Array{Complex{T}, 2}`: the STFT matrix, usually having (1+nfft/2) rows
* `samplerate::Real`: Sample rate of the resulting audio
* `windowsize::Int`: window size used for STFT. (default: nfft)
* `hopsize::Int`: number of frames between STFT columns (default: windowsize/4)
* `nfft::Int`: the FFT size (default: 2*(size(stft,1)-1))
* `window::Union{Function, AbstractVector, Nothing}`: function or vector that has the window used for STFT (default: uniform window)
"""
function istft(stft::Array{Complex{T}, 2},
                                   samplerate::Real,
                                   windowsize::Int = 2 * (size(stft, 1) - 1),
                                   hopsize::Int = windowsize >> 2;
                                   nfft::Int = windowsize,
                                   window::Union{Function, AbstractVector, Nothing} = hanning) where {T <: AbstractFloat}

    if windowsize > nfft
      error("window size should be less than or equal to nfft")
    end

    window::Vector{T} = map(T, w(window, windowsize))

    columns = size(stft, 2)
    audio = zeros(T, windowsize + hopsize * (columns - 1))
    weights = zeros(T, windowsize + hopsize * (columns - 1))
    spectrum = zeros(T, nfft)

    nbins = size(stft, 1)
    base = Base.unsafe_convert(Ptr{Complex{T}}, stft)
    stride = nbins * sizeof(Complex{T})

    @inbounds for i in 1:columns
        irfft!(spectrum, base + stride * (i-1), nfft)

        left = (i-1) * hopsize
        for j in 1:windowsize
            audio[left + j] += window[j] * spectrum[j]
            weights[left + j] += window[j]
        end
    end

    @inbounds for i in 1:length(audio)
        if weights[i] > eps(T)
            audio[i] /= weights[i]
        end
    end

    SampleBuf{T, 1}(audio, samplerate)
end


function istft(stft::Array{Complex{T}, 2},
                                   samplerate::Hertz,
                                   windowsize::Int = 2 * (size(stft, 1) - 1),
                                   hopsize::Int = windowsize >> 2;
                                   nfft::Int = windowsize,
                                   window::Union{Function, AbstractVector, Nothing} = hanning) where {T <: AbstractFloat}
    istft(stft,
            ustrip(samplerate),
            windowsize,
            hopsize,
            nfft,
            window)
end

"""
    istft(stft, samplerate, windowsize, hopsize; window)

Reconstruct a multichannel audio from its STFT.

# Arguments
* `stft::Array{Complex{T}, 3}`: the STFT matrix, usually having (1+nfft/2) rows
* `samplerate::Real`: Sample rate of the resulting audio
* `windowsize::Int`: window size used for STFT. (default: nfft)
* `hopsize::Int`: number of frames between STFT columns (default: windowsize/4)
* `nfft::Int`: the FFT size (default: 2*(size(stft,1)-1))
* `window::Union{Function, AbstractVector, Nothing}`: function or vector that has the window used for STFT (default: uniform window)
"""
function istft(stft::Array{Complex{T}, 3},
                                   samplerate::Real,
                                   windowsize::Int = 2 * (size(stft, 1) - 1),
                                   hopsize::Int = windowsize >> 2; kwargs...) where {T <: AbstractFloat}
    nchannels = size(stft, 3)
    buffers = Vector{SampleBuf}(undef, nchannels)

    # run ISTFT for each channel
    for i = 1:nchannels
        buffers[i] = istft(stft[:, :, i], samplerate, windowsize, hopsize; kwargs...)
    end

    nsamples = size(buffers[1], 1)
    audio = Array{T}(undef, nsamples, nchannels)

    for i = 1:nchannels
        audio[:, i] = buffers[i].data
    end

    # concatenate channels
    SampleBuf{T, 2}(
        audio,
        samplerate
    )
end


function istft(stft::Array{Complex{T}, 3},
                                   samplerate::Hertz,
                                   windowsize::Int = 2 * (size(stft, 1) - 1),
                                   hopsize::Int = windowsize >> 2; kwargs...) where {T <: AbstractFloat}
    istft(stft,
            ustrip(samplerate),
            windowsize,
            hopsize;
            kwargs...)
end

"""
    phase_vocoder(stft, speed, hopsize)

Phase vocoder. Given an STFT matrix, speed it up by a factor.
"""
function phase_vocoder(stft::Array{Complex{T}, 2},
                                           speed::Real,
                                           hopsize::Int = (size(stft, 1) - 1) >> 1) where {T <: AbstractFloat}

    nbins = size(stft, 1)
    nframes = size(stft, 2)
    nfft = (nbins - 1) * 2
    timesteps = 1f0:speed:nframes
    stretched = zeros(Complex{T}, nbins, length(timesteps))
    phase_advance = collect(range(0f0, T(π * hopsize), length=nbins))
    phase_acc = map(angle, stft[:,1])
    cis_phase = Array{Complex{T}}(undef,size(phase_acc))
    angle1 = Array{T}(undef, nbins)
    angle2 = Array{T}(undef, nbins)
    dphase = Array{T}(undef, nbins)
    mag = Array{T}(undef, nbins)
    twopi = T(2π)

    @inbounds for (t, step) in enumerate(timesteps)
        left = floor(Int, step)
        right = left + 1

        # weighting for linear magnitude interpolation
        alpha = step % 1f0
        if left < nframes
            for i = 1:nbins
                mag[i] = (1f0 - alpha) * abs(stft[i, left]) + alpha * abs(stft[i, right])
            end
        else
            for i = 1:nbins
                mag[i] = abs(stft[i, left])
            end
        end

        # store to output array
        cis!(cis_phase, phase_acc)
        cis_phase .= cis_phase .* mag
        for i in 1:nbins
            stretched[i, t] = cis_phase[i]
        end

        if t == size(stretched, 2)
            # nothing to do more
            break
        end

        for i in 1:nbins
            angle1[i] = map(angle, stft[i, left])
            angle2[i] = map(angle,stft[i, right])
            # compute phase advance
            dphase[i] = angle2[i] - angle1[i] - phase_advance[i]
            # wrap to -pi:pi range
            dphase[i] -= twopi * round(dphase[i] / twopi)
            # accumulate phase
            phase_acc[i] += phase_advance[i] + dphase[i]
        end
    end

    stretched
end

"""
    phase_vocoder(stft, speed, hopsize)

Phase vocoder. Given an STFT matrix, speed it up by a factor.
"""
function phase_vocoder(stft::Array{Complex{T}, 3},
                                           speed::Real,
                                           hopsize::Int = (size(stft, 1) - 1) >> 1) where {T <: AbstractFloat}
    mapslices(stft, dims=1:2) do stft
        phase_vocoder(stft, speed, hopsize)
    end
end

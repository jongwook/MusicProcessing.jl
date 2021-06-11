import FFTW.Frequencies

"""The Mel spectrogram"""
struct MelSpectrogram{T, F <: Union{Frequencies,AbstractRange}} <: DSP.Periodograms.TFR{T}
    power::Matrix{T}
    mels::F
    time::StepRangeLen{Float64}
end
heatmap(tfr::MelSpectrogram) = log10(tfr.power)
DSP.freq(tfr::MelSpectrogram) = tfr.mels
Base.time(tfr::MelSpectrogram) = tfr.time

"""A TFR subtype where each column represents an MFCC vector"""
struct MFCC{T, F <: Union{Frequencies,AbstractRange}} <: DSP.Periodograms.TFR{T}
    mfcc::Matrix{T}
    number::F
    time::StepRangeLen{Float64}
end
heatmap(tfr::MFCC) = tfr.mfcc
DSP.freq(tfr::MFCC) = tfr.number
Base.time(tfr::MFCC) = tfr.time

""""""
function hz_to_mel(frequencies::Union{F, AbstractArray{F}}) where {F <: Real}
    f_min = 0f0
    f_sp = 200f0 / 3

    mels = collect((frequencies .- f_min) ./ f_sp)

    min_log_hz = 1000f0
    min_log_mel = (min_log_hz - f_min) / f_sp
    logstep = log(6.4f0) / 27f0

    @inbounds for i = 1:length(mels)
        if frequencies[i] >= min_log_hz
            mels[i] = min_log_mel + log(frequencies[i] / min_log_hz) / logstep
        end
    end

    mels
end

""""""
function mel_to_hz(mels::Union{F, AbstractArray{F}}) where {F <: Real}
    f_min = 0f0
    f_sp = 200f0 / 3
    frequencies = collect(f_min .+ f_sp .* mels)

    min_log_hz = 1000f0
    min_log_mel = (min_log_hz - f_min) / f_sp
    logstep = log(6.4f0) / 27f0

    @inbounds for i = 1:length(frequencies)
        if mels[i] >= min_log_mel
            frequencies[i] = min_log_hz * exp(logstep * (mels[i] - min_log_mel))
        end
    end

    frequencies
end

""""""
function mel_frequencies(nmels::Int = 128, fmin::Real = 0.0f0, fmax::Real = 11025f0)
    min_mel = hz_to_mel(fmin)[1]
    max_mel = hz_to_mel(fmax)[1]

    mels = range(min_mel, max_mel, length=nmels)
    mel_to_hz(mels)
end

""""""
function mel(samplerate::Real, nfft::Int, nmels::Int = 128, fmin::Real = 0f0, fmax::Real = samplerate/2f0)
    second_dim = ((nfft >> 1) + 1)
    weights = zeros(Float32, nmels, second_dim)
    fftfreqs = fft_frequencies(samplerate, nfft)
    melfreqs = mel_frequencies(nmels + 2, fmin, fmax)
    enorm = 2f0 ./ (melfreqs[3:end] .- melfreqs[1:nmels])

    for i in 1:nmels
        lower = (fftfreqs .- melfreqs[i]) ./ (melfreqs[i+1] .- melfreqs[i])
        upper = (melfreqs[i+2] .- fftfreqs) ./ (melfreqs[i+2] .- melfreqs[i+1])

        for j in 1:second_dim
            weights[i, j] = max(0, min(lower[j], upper[j])) * enorm[i]
        end

    end

    weights
end

function melspectrogram(audio::SampleBuf{T, 1}, windowsize::Int = 1024, hopsize::Int = windowsize >> 2;
                        nmels::Int = 128, fmin::Real = 0f0, fmax::Real = audio.samplerate / 2f0) where T
    samplerate = audio.samplerate
    nfft = DSP.nextfastfft(windowsize)
    S = spectrogram(audio, windowsize, hopsize).power
    data = mel(samplerate, nfft, nmels, fmin, fmax) * S
    nframes = size(data, 2)
    MelSpectrogram(data, range(hz_to_mel(fmin)[1], hz_to_mel(fmax)[1], length=nmels), (0.0:nframes-1) * hopsize / samplerate)
end

function mfcc(audio::SampleBuf{T, 1}, windowsize::Int = 1024, hopsize::Int = windowsize >> 2;
              nmfcc::Int = 20, nmels::Int = 128, fmin::Real = 0f0, fmax::Real = audio.samplerate / 2f0) where T
    if nmfcc >= nmels
        error("number of mfcc components should be less than the number of mel frequency bins")
    end

    M = melspectrogram(audio, windowsize, hopsize; nmels = nmels, fmin = fmin, fmax = fmax)
    mfcc = dct(nmfcc, nmels) * power(M)

    for frame in 1:size(mfcc, 2)
        mfcc[:, frame] /= norm(mfcc[:, frame])
    end

    MFCC(mfcc, 1:nmfcc, time(M))
end

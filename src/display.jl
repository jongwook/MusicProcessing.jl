using .GLMakie
export waveplot, specplot, test

"""
    test()
"""

function test()
    x = range(0, 10, length=100)
    y = sin.(x)
    lines(x, y)
end

"""
    specplot(audio::Sample{T,N}, fs = 44100Hz)

Draws spectrogram of an audio

# Example

```julia
using GLMakie
audio_one_channel = SampleBuf(rand(1000), 10)
specplot(audio_one_channel,22100)
```
"""
function specplot(audio::SampleBuf{T,N}, fs = 44100) where {T,N}
    n = length(audio.data)
    nw = n÷50
    spec = spectrogram(mono(audio).data, nw, nw÷2; fs=fs)
    Makie.heatmap(spec.time, spec.freq, pow2db.(spec.power))
end
"""
    waveplot(audio::Sample{T,N}, fs = 44100Hz)

Draws waveplot of a audio

# Example

```julia
using GLMakie
audio_one_channel = SampleBuf(rand(1000), 10)
waveplot(audio_one_channel,48000)
```
"""
function waveplot(audio::SampleBuf{T,N}, fs = 44100) where {T,N}
    lines(0:1/fs:(length(mono(audio))-1)/fs, mono(audio))
end
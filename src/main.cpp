#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <fstream>
#include <numeric>

// Bibliotheken einbinden
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>
#include <structs.hpp>
#include <AudioFile.h>
#include <cxxopts.hpp>
using namespace std;

[[nodiscard]] inline bool ffmpegExists() {
    // Prüft ob ffmpeg im PATH ist
    return std::system("ffmpeg -version > /dev/null 2>&1") == 0;
}

inline bool convertWithFFmpeg(const std::string &wavPath, const std::string &targetPath) {
    const std::string cmd = "ffmpeg -y -i \"" + wavPath + "\" \"" + targetPath + "\"";
    const int ret = std::system((cmd + " 2>&1").c_str());
    return ret == 0;
}

// --- Main ---
int main(int argc, char *argv[]) {
    cxxopts::Options options("ImageToSound", "Convert image to sound");
    options.add_options("General options")
            ("h,help", "Print help")
            ("i,input", "Input image path", cxxopts::value<std::string>()->default_value("Silly_Cat_Character.jpg"))
            ("o,output", "Output sound path", cxxopts::value<std::string>()->default_value(""));

    options.add_options("Audio options")
            ("min-freq", "Minimum frequency", cxxopts::value<float>()->default_value("400.0"))
            ("max-freq", "Maximum frequency", cxxopts::value<float>()->default_value("14000.0"))
            ("samplerate", "Sample rate", cxxopts::value<int>()->default_value("44100"))
            ("duration-per-column", "Duration per image column in seconds",
             cxxopts::value<float>()->default_value("0.01"));

    options.add_options("Output options")
        ("f,format", "Convert WAV to another format using ffmpeg (e.g., mp3, flac)", cxxopts::value<std::string>()->default_value("wav"))
        ("keep-wav", "Keep intermediate WAV file when converting to another format");

    const cxxopts::ParseResult parsed = options.parse(argc, argv);
    if (parsed.count("help")) {
        std::cout << options.help() << std::endl;
        std::cout << "When using -f to convert to another format, except wav, ffmpeg must be installed and available in PATH." << std::endl;
        return 0;
    }

    // --- Pfade ---
    const string inputImagePath = parsed["input"].as<std::string>();
    if (inputImagePath.empty()) throw std::runtime_error("Input path cannot be empty");

    string outputSoundPath;
    if (parsed["output"].as<std::string>().empty()) {
        outputSoundPath = inputImagePath.substr(0, inputImagePath.find_last_of('.')) + ".wav";
    } else {
        outputSoundPath = parsed["output"].as<std::string>();
    }

    // Pfade vorbereiten
    const std::string outputFormat = parsed["format"].as<std::string>(); // z.B. "mp3", "flac", "wav"
    std::filesystem::path wavPath(outputSoundPath);
    std::filesystem::path finalOutputPath = wavPath; // default WAV

    // Falls ein anderes Format gewünscht wird, nur Dateiendung ändern
    if (outputFormat != "wav") {
        finalOutputPath.replace_extension(outputFormat);
    }

    bool keepWav = parsed.count("keep-wav") > 0;

    // --- Image laden ---
    Image img;
    if (!img.load(inputImagePath)) {
        cerr << "Failed to load image: " << inputImagePath << endl;
        return 1;
    }

    // --- Audio-Parameter ---
    AudioParams params;
    params.samplerate = parsed["samplerate"].as<int>();
    params.minFreq = parsed["min-freq"].as<float>();
    params.maxFreq = parsed["max-freq"].as<float>();
    params.durationPerColumn = parsed["duration-per-column"].as<float>();

    int samplesPerColumn = static_cast<int>(static_cast<float>(params.samplerate) * params.durationPerColumn);

    // --- Audio bauen ---
    int totalSamples = 0;
    vector<int> frameOffsets;
    for (auto &frame: img.frames) {
        int colSamples = max(256, samplesPerColumn);
        frameOffsets.push_back(totalSamples);
        totalSamples += frame.width * colSamples;
    }

    // Reserve finalAudio
    vector<float> finalAudio(totalSamples, 0.0f);
    for (size_t f = 0; f < img.frames.size(); ++f) {
        auto &frame = img.frames[f];
        int offset = frameOffsets[f];
        vector<float> frequencies(frame.height);
        vector<float> phaseInc(frame.height);
        for (int y = 0; y < frame.height; ++y) {
            frequencies[y] = params.minFreq * pow(params.maxFreq / params.minFreq,
                                                  static_cast<float>(frame.height - 1 - y) / static_cast<float>(
                                                      frame.height - 1));
            phaseInc[y] = M_PI * 2.0f * frequencies[y] / params.samplerate;
        }
        // Jede Spalte bekommt mindestens 256 Samples für scharfes Spektrogramm
        int colSamples = max(256, samplesPerColumn);

        // Hann-Window vorberechnen (reduziert Spectral Leakage)
        vector<float> hann(colSamples);
        for (int i = 0; i < colSamples; ++i) {
            hann[i] = 0.5f * static_cast<float>(1.0f - cos(2.0f * M_PI * i / (colSamples - 1)));
        }

        // Spalten parallel berechnen
        #pragma omp parallel for schedule(static) default(none) \
            shared(finalAudio, frame, frequencies, hann, offset, colSamples,params, phaseInc)
        for (int x = 0; x < frame.width; ++x) {
            vector<float> phase(frame.height, 0.0f); // Phase pro Frequenz
            vector<float> amp(frame.height); // Amplitude pro Frequenz
            // Amplitude vorberechnen
            for (int y = 0; y < frame.height; ++y) {
                float c = frame.pixels[y * frame.width + x] - 0.5f;
                amp[y] = c * fabs(c);
            }
            // Samples der Spalte berechnen
            for (int i = 0; i < colSamples; ++i) {
                float sample = 0.0f;
                for (int y = 0; y < frame.height; ++y) {
                    sample += amp[y] * sinf(phase[y]);
                    phase[y] += phaseInc[y];
                }
                sample *= hann[i];
                finalAudio[offset + x * colSamples + i] += sample;
            }
        }
    }

    float mean = std::accumulate(finalAudio.begin(), finalAudio.end(), 0.0f)
                 / static_cast<float>(finalAudio.size());

    for (auto &s: finalAudio)
        s -= mean;

    // Normalisieren
    float maxAmp = 0.0f;
    for (const float v: finalAudio) maxAmp = std::max(maxAmp, std::abs(v));
    if (maxAmp > 0.0f) {
        for (float &v: finalAudio) v /= maxAmp;
    }

    // --- WAV speichern ---
    AudioFile<float> wav;
    wav.setNumChannels(1);
    wav.setSampleRate(params.samplerate);
    wav.setNumSamplesPerChannel(static_cast<int>(finalAudio.size()));
    for (int i = 0; i < finalAudio.size(); ++i)
        wav.samples[0][i] = finalAudio[i];

    // Falls ein anderes Format gewünscht ist, erst als WAV speichern und dann konvertieren
    wav.save(outputSoundPath);
    cout << "Audio saved to " << outputSoundPath << endl;

    if (outputFormat != "wav") {
        if (!ffmpegExists()) {
            cerr << "ffmpeg is not installed or not found in PATH. Cannot convert to " << outputFormat << endl;
            return 1;
        }
        if (!convertWithFFmpeg(outputSoundPath, finalOutputPath.string())) {
            cerr << "Failed to convert WAV to " << outputFormat << endl;
            cerr << "Make sure ffmpeg supports the format." << endl;
            return 1;
        }
        // Original WAV löschen
        if (!keepWav) {
            std::remove(outputSoundPath.c_str());
        }
        cout << "Converted audio saved to " << finalOutputPath.string() << endl;
    }
    return 0;
}

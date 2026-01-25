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
#include <AudioFile.h>
#include <cxxopts.hpp>

// Für meine eigenen sachen
#include <structs.hpp>
#include <ungrouped_utility_funktions.hpp>
using namespace std;

// --- Funktionen ---

/**
 * @brief Berechnet die Frequenz für eine gegebene y-Position im Bild
 * @param y die y-Position (0 = unten, height-1 = oben)
 * @param height die Höhe des Bildes
 * @param params die Audio-Parameter
 * @return die Frequenz in Hz
 * @author Lupo
 */
float freqForY(const int y, const int height, const AudioParams &params) {
    if (height <= 1) // Vermeidung Division durch null
        return params.minFreq;

    const float t = static_cast<float>(height - 1 - y) /
                    static_cast<float>(height - 1);

    switch (params.scaleType) {
        case AudioParams::ScaleType::LINEAR:
            return params.minFreq +
                   t * (params.maxFreq - params.minFreq);

        case AudioParams::ScaleType::LOGARITHMIC:
            return params.minFreq *
                   std::pow(params.maxFreq / params.minFreq, t);

        case AudioParams::ScaleType::MEL: {
            const float melMin = hzToMel(params.minFreq);
            const float melMax = hzToMel(params.maxFreq);
            const float mel = melMin + t * (melMax - melMin);
            return melToHz(mel);
        }

        case AudioParams::ScaleType::BARK: {
            const float barkMin = hzToBark(params.minFreq);
            const float barkMax = hzToBark(params.maxFreq);
            const float bark = barkMin + t * (barkMax - barkMin);
            return barkToHz(bark);
        }
    }
    return params.minFreq; // unreachable, aber Compiler happy
}

[[nodiscard]] bool save_wav_file(const bool useStereo, const std::string &outputSoundPath,
    const std::vector<float> &finalAudio, const std::vector<float> &finalAudioL, const std::vector<float> &finalAudioR) {
    try {
        // --- WAV speichern ---
        AudioFile<float> wav;
        if (useStereo) {
            wav.setNumChannels(2);
            wav.setNumSamplesPerChannel(static_cast<int>(finalAudioL.size()));
            for (size_t i = 0; i < finalAudioL.size(); ++i) {
                wav.samples[0][i] = finalAudioL[i];
                wav.samples[1][i] = finalAudioR[i];
            }
        } else {
            wav.setNumChannels(1);
            wav.setNumSamplesPerChannel(static_cast<int>(finalAudio.size()));
            for (size_t i = 0; i < finalAudio.size(); ++i)
                wav.samples[0][i] = finalAudio[i];
        }


        // Falls ein anderes Format gewünscht ist, erst als WAV speichern und dann konvertieren
        wav.save(outputSoundPath);
        cout << "Audio saved to " << outputSoundPath << endl;
        return true;
    } catch (std::exception &e) {
        cerr << "Error saving WAV file: " << e.what() << endl;
        return false;
    }
}

[[nodiscard]] bool convert_file(const std::string &outputFormat, const std::string &outputSoundPath, const bool keepWav,
                                              const std::filesystem::path &finalOutputPath = "./" /*Dummy Parameter*/) {
    clog << "Converting WAV to " << outputFormat << " using ffmpeg..." << endl;
    if (!ffmpegExists()) {
        cerr << "ffmpeg is not installed or not found in PATH. Cannot convert to " << outputFormat << endl;
        cerr << "Will keep the WAV file at " << outputSoundPath << endl;
        return false;
    }
    if (!convertWithFFmpeg(outputSoundPath, finalOutputPath.string())) {
        cerr << "Failed to convert WAV to " << outputFormat << endl;
        cerr << "Make sure ffmpeg supports the format." << endl;
        cerr << "WAV file is kept at " << outputSoundPath << endl;
        cerr << "See ffmpeg.log for details." << endl;
        return false;
    }
    // Original WAV löschen
    if (!keepWav) {
        std::remove(outputSoundPath.c_str());
        cout << "Intermediate WAV file deleted." << endl;
    }
    cout << "Converted audio saved to " << finalOutputPath.string() << endl;
    return true;
}

enum class StereoNorm {
    LINKED, // beide Kanäle gemeinsam (empfohlen)
    INDEPENDENT // L und R getrennt
};


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
            ("duration-per-column", "Duration per image column in seconds", cxxopts::value<float>()->default_value("0.01"))
            ("scale", "The scale type", cxxopts::value<std::string>()->default_value("logarithmic"))
            ("gamma", "Gamma correction", cxxopts::value<float>()->default_value("1.0"))
            ("norm", "Normalization: peak|rms", cxxopts::value<std::string>()->default_value("peak"))
            ("stereo", "Enable stereo output (default is mono)");

    options.add_options("Output options")
            ("f,format", "Convert WAV to another format using ffmpeg (e.g., mp3, flac)",
             cxxopts::value<std::string>()->default_value("wav"))
            ("keep-wav", "Keep intermediate WAV file when converting to another format");

    const cxxopts::ParseResult parsed = options.parse(argc, argv);
    if (parsed.count("help")) {
        std::cout << options.help() << std::endl;
        std::cout << "When using -f to convert to another format, "
                  << " except wav, ffmpeg must be installed and available in PATH."
                  << std::endl;
        return 0;
    }

    // --- Pfade ---
    const string inputImagePath = parsed["input"].as<std::string>(); // Path to input image
    if (inputImagePath.empty()) throw std::runtime_error("Input path cannot be empty");

    string outputSoundPath; // Path to output sound
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
    Image img; // Das bild das in audio umgewandelt werden soll
    if (!img.load(inputImagePath)) {
        cerr << "Failed to load image: " << inputImagePath << endl;
        return 1;
    }

    // --- Audio-Parameter ---
    AudioParams params; // Die audioparameter für die sound erstellung
    params.samplerate = parsed["samplerate"].as<int>();
    params.minFreq = parsed["min-freq"].as<float>();
    params.maxFreq = parsed["max-freq"].as<float>();
    params.durationPerColumn = parsed["duration-per-column"].as<float>();
    params.gamma = parsed["gamma"].as<float>();

    auto toLower = [](std::string s) -> std::string {
        ranges::transform(s, s.begin(), ::tolower);
        return s;
    };

    // ReSharper disable once CppTooWideScopeInitStatement
    const std::string scaleStr = toLower(parsed["scale"].as<std::string>());
    if (scaleStr == "linear") {
        params.scaleType = AudioParams::ScaleType::LINEAR;
    } else if (scaleStr == "logarithmic" || scaleStr == "log") {
        params.scaleType = AudioParams::ScaleType::LOGARITHMIC;
    } else if (scaleStr == "mel") {
        params.scaleType = AudioParams::ScaleType::MEL;
    } else if (scaleStr == "bark") {
        params.scaleType = AudioParams::ScaleType::BARK;
    } else {
        throw std::runtime_error("Invalid scale type: " + scaleStr + " (use linear|log|mel|bark)");
    }

    std::string normStr = toLower(parsed["norm"].as<std::string>()); // User input string for normalization
    enum class NormType { PEAK, RMS } norm; // Normalization type to use

    if (normStr == "peak") {
        norm = NormType::PEAK;
    } else if (normStr == "rms") {
        norm = NormType::RMS;
    } else {
        throw std::runtime_error("Invalid normalization type: " + normStr + " (use peak|rms)");
    }


    int samplesPerColumn = static_cast<int>(static_cast<float>(params.samplerate) * params.durationPerColumn); // Samples pro Spalte

    // --- Audio bauen ---
    int totalSamples = 0; // Gesamtanzahl der Samples
    vector<int> frameOffsets; // Offset für jeden Frame im finalen Audio
    frameOffsets.reserve(img.frames.size()); // Speicher reservieren für Frame-Offsets
    for (auto &frame: img.frames) {
        int colSamples = max(256, samplesPerColumn); // Jede Spalte bekommt mindestens 256 Samples für scharfes Spektrogramm
        frameOffsets.push_back(totalSamples);            // Offset für diesen Frame
        totalSamples += frame.width * colSamples;        // Samples für diesen Frame hinzufügen
    }

    clog << "Generating audio: " << totalSamples << " samples at " << params.samplerate << " Hz" << endl;
    clog << "This may take a while depending on image/gif size and number of frames..." << endl;

    // Reserve finalAudio
    bool useStereo = parsed.count("stereo") > 0; // Ob Stereo verwendet werden soll oder Mono

    vector<float> finalAudio;               // Für Mono
    vector<float> finalAudioL, finalAudioR; // Für Stereo

    if (useStereo) {
        finalAudioL.resize(totalSamples, 0.0f);
        finalAudioR.resize(totalSamples, 0.0f);
    } else {
        finalAudio.resize(totalSamples, 0.0f);
    }

    // Jeden Frame durchgehen
    for (size_t f = 0; f < img.frames.size(); ++f) {
        auto &frame = img.frames[f];            // aktueller Frame
        int offset = frameOffsets[f];           // Frequenzen und Phaseninkremente für jede Zeile vorberechnen
        vector<float> frequencies(frame.height);// Frequenzen für jede Zeile
        vector<float> phaseInc(frame.height);   // Phaseninkrement pro Sample für jede Zeile

        // Frequenzen und Phaseninkremente berechnen
        for (int y = 0; y < frame.height; ++y) {
            frequencies[y] = freqForY(y, frame.height, params); // Frequenz für diese Zeile
            phaseInc[y] = M_PI * 2.0f * frequencies[y] / params.samplerate; // Phaseninkrement pro Sample
        }

        // Jede Spalte bekommt mindestens 256 Samples für scharfes Spektrogramm
        int colSamples = max(256, samplesPerColumn);

        // Hann-Window vorberechnen (reduziert Spectral Leakage)
        vector<float> hann(colSamples);
        for (int i = 0; i < colSamples; ++i) {
            hann[i] = 0.5f * static_cast<float>(1.0f - cos(2.0f * M_PI * i / (colSamples - 1)));
        }

        // Spalten parallel berechnen
        #pragma omp parallel for schedule(static) default(none) shared(finalAudio, finalAudioL, finalAudioR, frame,\
            frequencies, hann, offset, colSamples, params, phaseInc, useStereo)
        // Jede Spalte durchgehen
        for (int x = 0; x < frame.width; ++x) {
            // Panning berechnen
            float pan = (frame.width > 1)
                            ? 2.0f * (static_cast<float>(x) / static_cast<float>(frame.width - 1)) - 1.0f
                            : 0.0f;

            // Equal Power Panning
            float panL = 1.0f, panR = 1.0f; // Default für Mono
            if (useStereo) {
                panL = std::sqrt(0.5f * (1.0f - pan));
                panR = std::sqrt(0.5f * (1.0f + pan));
            }

            vector<float> phase(frame.height);
            vector<float> amp(frame.height);

            // Amplituden und Startphasen für jede Zeile berechnen
            for (int y = 0; y < frame.height; ++y) {
                float c = frame.pixels[y * frame.width + x] - 0.5f;
                float sign = (c >= 0.0f) ? 1.0f : -1.0f;
                amp[y] = sign * std::pow(std::abs(c), params.gamma);
                phase[y] = phaseInc[y] * static_cast<float>(x) * static_cast<float>(colSamples); // optional, aber gut
            }

            // Samples für diese Spalte generieren
            for (int i = 0; i < colSamples; ++i) {
                float sample = 0.0f; // Sample für diese Position
                for (int y = 0; y < frame.height; ++y) { // Alle Zeilen addieren
                    sample += amp[y] * sinf(phase[y]);
                    phase[y] += phaseInc[y];
                }
                sample *= hann[i]; // Hann-Window anwenden

                int idx = offset + x * colSamples + i; // Index im finalen Audio
                if (useStereo) { // Stereo
                    finalAudioL[idx] += sample * panL;
                    finalAudioR[idx] += sample * panR;
                } else { // Mono
                    finalAudio[idx] += sample; // Mono
                }
            }
        }
    }
    clog << "Audio generation completed." << endl;
    clog << "Post-processing audio..." << endl;

    // --- DC entfernen ---

    // DC-Offset entfernen
    auto removeDC = [](vector<float> &buf) -> void {
        const float mean = std::accumulate(buf.begin(), buf.end(), 0.0f) / static_cast<float>(buf.size());
        for (auto &s: buf) s -= mean;
    };

    if (useStereo) {
        removeDC(finalAudioL);
        removeDC(finalAudioR);
    } else {
        removeDC(finalAudio);
    }

    // --- Normalisieren ---
    if (useStereo) { // Stereo
        if (norm == NormType::RMS) { // RMS
            double sumSq = 0.0; // Summe der Quadrate
            for (size_t i = 0; i < finalAudioL.size(); ++i) { // Beide Kanäle zusammen betrachten (wir müssen nur die größe von einem Kanal kennen, sind ja gleich groß)
                sumSq += finalAudioL[i] * finalAudioL[i];
                sumSq += finalAudioR[i] * finalAudioR[i];
            }
            auto rms = static_cast<float>(std::sqrt(sumSq / (2.0 * static_cast<float>(finalAudioL.size())))); // RMS über beide Kanäle
            float gain = 0.1f / std::max(rms, 1e-6f);   // Verstärkung berechnen
            for (size_t i = 0; i < finalAudioL.size(); ++i) { // Beide Kanäle normalisieren
                finalAudioL[i] *= gain;
                finalAudioR[i] *= gain;
            }
        } else {
            // PEAK
            float maxAmp = 0.0f; // Maximale Amplitude über beide Kanäle finden
            for (size_t i = 0; i < finalAudioL.size(); ++i) // Beide Kanäle durchgehen und Maximum finden
                maxAmp = std::max({maxAmp, std::abs(finalAudioL[i]), std::abs(finalAudioR[i])});
            if (maxAmp > 0.0f) { // Normalisieren, falls maxAmp > 0
                float inv = 1.0f / maxAmp; // Invers der maximalen Amplitude
                for (size_t i = 0; i < finalAudioL.size(); ++i) { // Beide Kanäle normalisieren
                    finalAudioL[i] *= inv;
                    finalAudioR[i] *= inv;
                }
            }
        }
    } else {
        if (norm == NormType::RMS) {
            double sumSq = 0.0;
            for (float s: finalAudio) sumSq += s * s;
            auto rms = static_cast<float>(std::sqrt(sumSq / static_cast<float>(finalAudio.size()))); // RMS berechnen
            float gain = 0.1f / std::max(rms, 1e-6f); // Verstärkung berechnen
            for (float &s: finalAudio) s *= gain; // Normalisieren
        } else {
            // PEAK
            float maxAmp = *ranges::max_element(finalAudio, [](const float a, const float b) -> bool { // Maximum der absoluten Werte finden
                return std::abs(a) < std::abs(b);
            });
            if (maxAmp > 0.0f) { // Normalisieren, falls maxAmp > 0
                float inv = 1.0f / maxAmp;
                for (auto &s: finalAudio) s *= inv;
            }
        }
    }


    clog << "Audio post-processing completed." << endl;
    clog << "Saving audio to file..." << endl;

    // --- WAV speichern ---
    if (!save_wav_file(useStereo, wavPath.string(), finalAudio, finalAudioL, finalAudioR)) { // WAV speichern
        cerr << "Could not save WAV file to " << wavPath.string() << endl;
        return 1;
    }

    // --- Falls gewünscht, in anderes Format konvertieren ---
    if (outputFormat != "wav") {
        if (!convert_file(outputFormat, wavPath.string(), keepWav, finalOutputPath)) { // Konvertieren
            cerr << "Could not convert file " << wavPath.string() << endl;
            cerr << "Wav file is kept." << endl;
            return 1;
        }
    }

    clog << "Good bye! :3" << endl;
    return 0;
}

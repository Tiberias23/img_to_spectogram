#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <numeric>
#include <random>

// Bibliotheken einbinden
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>
#include <AudioFile.h> // TODO: evtl. selber schreiben, da die lib viel bietet aber wir nur wav speichern brauchen
#include <cxxopts.hpp>

// Für meine eigenen sachen
#include <structs.hpp> // für Image, ImageFrame, AudioParams, NormType und StereoNorm
#include <ungrouped_utility_funktions.hpp> // für freqForY, convert_file, ffmpegExists, convertWithFFmpeg, hzToMel, melToHz, hzToBark, barkToHz
using namespace std;


/**
 * @brief Speichert das Audio als WAV-Datei
 * @param useStereo ob Stereo verwendet wird
 * @param outputSoundPath der Pfad zur Ausgabedatei
 * @param finalAudio das finale Audio (Mono)
 * @param finalAudioL das finale Audio links (Stereo)
 * @param finalAudioR das finale Audio rechts (Stereo)
 * @param samplerate die rate der samples
 * @return true bei Erfolg, false bei Fehler
 * @author Lupo
 */
[[nodiscard]] inline bool save_wav_file(const bool useStereo, const std::string &outputSoundPath,
    const std::vector<float> &finalAudio, const std::vector<float> &finalAudioL, const std::vector<float> &finalAudioR, const int samplerate) {
    try {
        // --- WAV speichern ---
        AudioFile<float> wav;
        wav.setSampleRate(samplerate);
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

        if (!wav.save(outputSoundPath, AudioFileFormat::Wave)) {
            throw std::runtime_error("Failed to save file");
        }
        std::cout << "Audio saved to " << outputSoundPath << std::endl;
        return true;
    } catch (std::exception &e) {
        std::cerr << "Error saving WAV file: " << e.what() << std::endl;
        return false;
    }
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

    // konvert string to lower
    auto toLower = [](std::string s) -> std::string {
        ranges::transform(s, s.begin(), ::tolower);
        return s;
    };

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
    const std::string outputFormat = toLower(parsed["format"].as<std::string>()); // z.B. "mp3", "flac", "wav"
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


    // Samples pro Spalte
    int samplesPerColumn = static_cast<int>(static_cast<float>(params.samplerate) * params.durationPerColumn);

    // --- Audio bauen ---
    int totalSamples = 0; // Gesamtanzahl der Samples
    vector<int> frameOffsets; // Offset für jeden Frame im finalen Audio
    frameOffsets.reserve(img.frames.size()); // Speicher reservieren für Frame-Offsets
    for (auto &frame: img.frames) {
        // Jede Spalte bekommt mindestens 256 Samples für scharfes Spektrogramm
        int colSamples = max(256, samplesPerColumn);
        frameOffsets.push_back(totalSamples);            // Offset für diesen Frame
        totalSamples += frame.width * colSamples;        // Samples für diesen Frame hinzufügen
    }

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

    clog << "Generating audio: " << totalSamples << " samples at " << params.samplerate << " Hz" << endl;
    clog << "This may take a while depending on image/gif size and number of frames..." << endl;

    constexpr int LUT_SIZE = 4096;
    static vector<float> gammaLUT;
    static float lastGamma = -1.0f;

    if (gammaLUT.empty() || lastGamma != params.gamma) {
        gammaLUT.resize(LUT_SIZE);
        for (int i = 0; i < LUT_SIZE; ++i) {
            float v = static_cast<float>(i) / (LUT_SIZE - 1);
            gammaLUT[i] = std::pow(v, params.gamma);
        }
        lastGamma = params.gamma;
    }


    // TODO: Mach den scheiß schneller
    // TODO: eventuell was gegen diesen mini teil bei jedem frame links unten machen
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

        vector<float> sinInc(frame.height);
        vector<float> cosInc(frame.height);
        for (int y = 0; y < frame.height; ++y) {
            sinInc[y] = sinf(phaseInc[y]);
            cosInc[y] = cosf(phaseInc[y]);
        }

        // Jede Spalte bekommt mindestens 256 Samples für scharfes Spektrogramm
        int colSamples = max(256, samplesPerColumn);

        // Hann-Window vorberechnen (reduziert Spectral Leakage)
        vector<float> hann(colSamples);
        for (int i = 0; i < colSamples; ++i) {
            hann[i] = 0.5f * static_cast<float>(1.0f - cos(2.0f * M_PI * i / (colSamples - 1)));
        }

        const bool useWindow = img.frames.size() > 1; // Check ob Windowing verwendet werden soll (bei GIFs ja, bei Einzelbildern nein)

        // Spalten parallel berechnen
#pragma omp parallel default(none) \
        shared(finalAudio, finalAudioL, finalAudioR, frame, frequencies, hann, offset, colSamples, params, phaseInc, useStereo, img, useWindow, sinInc, cosInc, gammaLUT)
        // ReSharper disable once CppDFAUnreadVariable
        {
            vector<float> amp(frame.height);
            vector<float> sinv(frame.height);
            vector<float> cosv(frame.height);

            // Jede Spalte durchgehen
            #pragma omp for schedule(static, 16)
            for (int x = 0; x < frame.width; ++x) {
                // Panning berechnen
                float pan;
                if ((frame.width > 1))
                    pan = 2.0f * (static_cast<float>(x) / static_cast<float>(frame.width - 1)) - 1.0f;
                else
                    pan = 0.0f;

                // Equal Power Panning
                float panL = 1.0f, panR = 1.0f; // Default für Mono
                if (useStereo) {
                    panL = std::sqrt(0.5f * (1.0f - pan));
                    panR = std::sqrt(0.5f * (1.0f + pan));
                }

                // Amplituden und Startphasen für jede Zeile berechnen
                for (int y = 0; y < frame.height; ++y) {
                    float c = frame.pixels[y * frame.width + x]; // Pixelwert an (x,y)
                    if (params.scaleType == AudioParams::ScaleType::LINEAR || img.frames.size() > 1 /* GIF erkennen: mehrere Frames*/) {
                        // alte Sign-Pow Logik für LINEAR und GIFs
                        c -= 0.5f;
                        float sign = (c >= 0.0f) ? 1.0f : -1.0f;
                        float ac = std::abs(c);
                        int li = std::min(static_cast<int>(ac * (LUT_SIZE - 1)), LUT_SIZE - 1);
                        amp[y] = sign * gammaLUT[li];
                    } else {
                        // LOG / MEL / BARK für JPGs
                        c = std::max(c, 1e-4f);
                        amp[y] = std::pow(c, params.gamma);
                    }

                    float startPhase = phaseInc[y] * static_cast<float>(x) * static_cast<float>(colSamples);
                    sinv[y] = sinf(startPhase);
                    cosv[y] = cosf(startPhase);

                }


                // Samples für diese Spalte generieren
                const int colSamplesFrame = colSamples; // Fix für diese Schleife
                for (int i = 0; i < colSamplesFrame; ++i) {
                    float sample = 0.0f; // Sample für diese Position
                    for (int y = 0; y < frame.height; ++y) { // Alle Zeilen addieren
                        sample += amp[y] * sinv[y];

                        float s = sinv[y];
                        float c = cosv[y];
                        sinv[y] = s * cosInc[y] + c * sinInc[y];
                        cosv[y] = c * cosInc[y] - s * sinInc[y];
                    }
                    if (useWindow) sample *= hann[i]; // Hann-Window anwenden

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
            // beide Kanäle zusammen betrachten (wir müssen nur die größe von einem Kanal kennen, sind ja gleich groß)
            for (size_t i = 0; i < finalAudioL.size(); ++i) {
                sumSq += finalAudioL[i] * finalAudioL[i];
                sumSq += finalAudioR[i] * finalAudioR[i];
            }
            // RMS über beide Kanäle
            auto rms = static_cast<float>(std::sqrt(sumSq / (2.0 * static_cast<float>(finalAudioL.size()))));
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
            // Maximum der absoluten Werte finden
            float maxAmp = *ranges::max_element(finalAudio, [](const float a, const float b) -> bool {
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
    // WAV speichern
    if (!save_wav_file(useStereo, wavPath.string(), finalAudio, finalAudioL, finalAudioR, params.samplerate)) {
        cerr << "Could not save WAV file to " << wavPath.string() << endl;
        return 1;
    }
    clog << "WAV file saved successfully." << endl;

    // --- Falls gewünscht, in anderes Format konvertieren ---
    if (outputFormat != "wav") {
        if (!convert_file(outputFormat, wavPath.string(), keepWav, finalOutputPath)) { // Konvertieren
            cerr << "Could not convert file " << wavPath.string() << endl;
            cerr << "Wav file is kept." << endl;
            return 1;
        }
        clog << "File conversion completed Successfully." << endl;
    }

    clog << "Good bye! :3" << endl;
    return 0;
}

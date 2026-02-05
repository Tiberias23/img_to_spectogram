#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <numeric>

// Include library's
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>
#include <AudioFile.h> // TODO: evtl. selber schreiben, da die lib viel bietet aber wir nur wav speichern brauchen
#include <cxxopts.hpp>

// Include the stuff form the other files
#include <structs.hpp> // for Image, ImageFrame, AudioParams, NormType und StereoNorm
#include <ungrouped_utility_funktions.hpp> // for freqForY, convert_file, ffmpegExists, convertWithFFmpeg, hzToMel, melToHz, hzToBark, barkToHz
using namespace std;


/**
 * @brief Save WAV file
 * @param useStereo whether to use stereo or mono
 * @param outputSoundPath the output path
 * @param finalAudio the finale Audio (Mono)
 * @param finalAudioL the finale Audio left (Stereo)
 * @param finalAudioR the finale Audio right (Stereo)
 * @param samplerate the samplerate
 * @return true on success, false on failure
 * @author Lupo
 */
[[nodiscard]] inline bool save_wav_file(const bool useStereo, const std::string &outputSoundPath,
    const std::vector<float> &finalAudio, const std::vector<float> &finalAudioL, const std::vector<float> &finalAudioR, const int samplerate) {
    try {
        // --- Save WAV ---
        AudioFile<float> wav;           // AudioFile object
        wav.setSampleRate(samplerate);  // Set sample rate
        if (useStereo) {
            wav.setNumChannels(2); // Stereo
            wav.setNumSamplesPerChannel(static_cast<int>(finalAudioL.size())); // Set number of samples
            for (size_t i = 0; i < finalAudioL.size(); ++i) { // Write samples to both channels
                wav.samples[0][i] = finalAudioL[i];
                wav.samples[1][i] = finalAudioR[i];
            }
        } else {
            wav.setNumChannels(1); // Mono
            wav.setNumSamplesPerChannel(static_cast<int>(finalAudio.size())); // Set number of samples
            for (size_t i = 0; i < finalAudio.size(); ++i) // Write samples to mono channel
                wav.samples[0][i] = finalAudio[i];
        }

        if (!wav.save(outputSoundPath, AudioFileFormat::Wave)) { // Save WAV file
            throw std::runtime_error("Failed to save file");     // Throw error on failure
        }
        std::cout << "Audio saved to " << outputSoundPath << std::endl;
        return true;
    } catch (std::exception &e) {
        std::cerr << "Error saving WAV file: " << e.what() << std::endl; // Print error message
        return false; // Return false on failure
    }
}


// --- Main ---
int main(int argc, char *argv[]) {
    cxxopts::Options options("ImageToSound", "Convert image to sound"); // Command line options
    options.add_options("General options")
            ("h,help", "Print help")
            ("i,input", "Input image path default: ./Silly_Cat_Character.jpg", cxxopts::value<std::string>()->default_value("Silly_Cat_Character.jpg"))
            ("o,output", "Output sound path default ./INPUT_FILENAME.wav", cxxopts::value<std::string>()->default_value(""));

    options.add_options("Audio options")
            ("min-freq", "Minimum frequency default 400.0", cxxopts::value<float>()->default_value("400.0"))
            ("max-freq", "Maximum frequency default 14000.0", cxxopts::value<float>()->default_value("14000.0"))
            ("samplerate", "Sample rate", cxxopts::value<int>()->default_value("44100"))
            ("duration-per-column", "Duration per image column in seconds", cxxopts::value<float>()->default_value("0.01"))
            ("scale", "The scale type (linear | logarithmic | mel | bark)", cxxopts::value<std::string>()->default_value("logarithmic"))
            ("gamma", "Gamma correction", cxxopts::value<float>()->default_value("1.0"))
            ("gama-size,gama-look-up-table-size", "Size of the gamma look-up table (only for linear scale)", cxxopts::value<int>()->default_value("4096"))
            ("norm", "Normalization: peak|rms", cxxopts::value<std::string>()->default_value("peak"))
            ("stereo", "Enable stereo output (default is mono)");

    options.add_options("Output options")
            ("f,format", "Convert WAV to another format using ffmpeg (e.g., mp3, flac)",
             cxxopts::value<std::string>()->default_value("wav"))
            ("keep-wav", "Keep intermediate WAV file when converting to another format");

    const cxxopts::ParseResult parsed = options.parse(argc, argv); // Parsed arguments
    if (parsed.count("help")) {
        std::cout << options.help() << std::endl;
        std::cout << "Note:\n\tWhen using -f to convert to another format, "
                  << " except wav, ffmpeg must be installed and available in PATH."
                  << std::endl;
        return 0;
    }

    // konvert a string to its lowercase representation
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
    const std::string outputFormat = toLower(parsed["format"].as<std::string>()); // the output format e.g. "mp3", "flac", "wav"
    std::filesystem::path wavPath(outputSoundPath);     // Path to wav file (temporary if converting)
    std::filesystem::path finalOutputPath = wavPath;    // the final output path

    // if the output format is not wav, change the final output path extension accordingly
    if (outputFormat != "wav") {
        finalOutputPath.replace_extension(outputFormat);
    }

    bool keepWav = parsed.count("keep-wav") > 0; // whether to keep the intermediate wav file

    // --- Image laden ---
    Image img; // The image object to hold the loaded image
    if (!img.load(inputImagePath)) {
        cerr << "Failed to load image: " << inputImagePath << endl;
        return 1;
    }

    // --- Audio-Parameter ---
    AudioParams params; // The audio parameters
    params.samplerate = parsed["samplerate"].as<int>(); // Sample rate
    params.minFreq = parsed["min-freq"].as<float>();    // Minimum frequency
    params.maxFreq = parsed["max-freq"].as<float>();    // Maximum frequency
    params.durationPerColumn = parsed["duration-per-column"].as<float>(); // Duration per image column
    params.gamma = parsed["gamma"].as<float>();         // Gamma correction

    const std::string scaleStr = toLower(parsed["scale"].as<std::string>()); // User input string for scale type
    if (scaleStr == "linear" || scaleStr == "lin") {
        params.scaleType = AudioParams::ScaleType::LINEAR;      // Linear scale
    } else if (scaleStr == "logarithmic" || scaleStr == "log") {
        params.scaleType = AudioParams::ScaleType::LOGARITHMIC; // Logarithmic scale
    } else if (scaleStr == "mel") {
        params.scaleType = AudioParams::ScaleType::MEL;         // Mel scale
    } else if (scaleStr == "bark") {
        params.scaleType = AudioParams::ScaleType::BARK;        // Bark scale
    } else {
        throw std::runtime_error("Invalid scale type: " + scaleStr + " (use linear|logarithmic|mel|bark)"); // Invalid scale type throw error
    }

    std::string normStr = toLower(parsed["norm"].as<std::string>()); // User input string for normalization
    enum class NormType { PEAK, RMS } norm; // Normalization type to use Normalisieren, falls maxAmp > 0
    if (normStr == "peak") {        // Peak normalization
        norm = NormType::PEAK;
    } else if (normStr == "rms") {  // RMS normalization
        norm = NormType::RMS;
    } else {
        throw std::runtime_error("Invalid normalization type: " + normStr + " (use peak|rms)"); // Invalid normalization type throw error
    }

    // Samples per column
    int samplesPerColumn = static_cast<int>(static_cast<float>(params.samplerate) * params.durationPerColumn);

    int lookupTableSize = parsed["gama-look-up-table-size"].as<int>(); // Size of the gamma look-up table
    if (lookupTableSize <= 0) {
        throw std::runtime_error("Gamma look-up table size must be positive");
    }

    // --- Audio bauen ---
    int totalSamples = 0;     // Total sample count
    vector<int> frameOffsets; // Offset for each frame in the final audio
    frameOffsets.reserve(img.frames.size()); // Reserve space for frame offsets
    for (auto &frame: img.frames) {
        // Every column gets at least 256 samples for a sharp spectrogram
        int colSamples = max(256, samplesPerColumn);
        frameOffsets.push_back(totalSamples);            // Offset for this frame
        totalSamples += frame.width * colSamples;        // add samples for this frame
    }

    // Reserve finalAudio
    bool useStereo = parsed.count("stereo") > 0; // Whether to use stereo output or mono output

    vector<float> finalAudio;               // For Mono
    vector<float> finalAudioL, finalAudioR; // For Stereo

    // Reserve space according to stereo/mono
    if (useStereo) {
        finalAudioL.resize(totalSamples, 0.0f);
        finalAudioR.resize(totalSamples, 0.0f);
    } else {
        finalAudio.resize(totalSamples, 0.0f);
    }

    clog << "Generating audio: " << totalSamples << " samples at " << params.samplerate << " Hz" << endl;
    clog << "This may take a while depending on image/gif size and number of frames..." << endl;

    const int LUT_SIZE = lookupTableSize;   // Size of the gamma lookup table
    static vector<float> gammaLUT;          // Gamma lookup table
    static float lastGamma = -1.0f;         // Last used gamma value (to check if we need to recalculate LUT)

    // Gamma-LUT vorberechnen (für LINEAR Skalierung)
    if (gammaLUT.empty() || lastGamma != params.gamma) {
        gammaLUT.resize(LUT_SIZE); // Resize LUT
        for (int i = 0; i < LUT_SIZE; ++i) { // Fill LUT
            float v = static_cast<float>(i) / static_cast<float>(LUT_SIZE - 1); // Normalized value [0,1]
            gammaLUT[i] = std::pow(v, params.gamma);      // Apply gamma correction
        }
        lastGamma = params.gamma; // Update last gamma
    }


    // TODO: eventuell was gegen diesen mini teil bei jedem frame links unten machen
    // Go through each frame
    for (size_t f = 0; f < img.frames.size(); ++f) {
        auto &frame = img.frames[f];            // Current frame
        int offset = frameOffsets[f];           // Offset in final audio for this frame
        vector<float> frequencies(frame.height);// Frequencies for each row
        vector<float> phaseInc(frame.height);   // Phase increments for each row

        // Calculate frequencies and phase increments for each row
        for (int y = 0; y < frame.height; ++y) {
            frequencies[y] = freqForY(y, frame.height, params);             // Frequency for this row
            phaseInc[y] = M_PI * 2.0f * frequencies[y] / params.samplerate; // Phase increment for this row
        }

        vector<float> sinInc(frame.height);         // Sine of phase increments
        vector<float> cosInc(frame.height);         // Cosine of phase increments
        for (int y = 0; y < frame.height; ++y) {    // Precompute sine and cosine of phase increments
            sinInc[y] = sinf(phaseInc[y]);          // Sine of phase increment
            cosInc[y] = cosf(phaseInc[y]);          // Cosine of phase increment
        }

        // Samples per column (at leste 256 for a sharp spectrogram)
        int colSamples = max(256, samplesPerColumn);

        vector<float> hann(colSamples);         // Hann window
        for (int i = 0; i < colSamples; ++i) {  // Precalculate Hann window for this column
            hann[i] = 0.5f * static_cast<float>(1.0f - cos(2.0f * M_PI * i / (colSamples - 1))); // Hann window formula
        }

        // Check whether to use Windowing (for GIFs yes, for single image's no)
        const bool useWindow = img.frames.size() > 1;

        // multithread the calculation with OpenMP
#pragma omp parallel default(none) shared(finalAudio, finalAudioL, finalAudioR, frame, frequencies, hann, offset,\
    colSamples, params, phaseInc, useStereo, img, useWindow, sinInc, cosInc, gammaLUT, LUT_SIZE)
        {
            vector<float> amp(frame.height);                // Amplitudes for each row
            vector<float> sine_values(frame.height);        // Sine values for each row
            vector<float> cosine_values(frame.height);      // Cosine values for each row

            // Go through each column
            #pragma omp for schedule(static, 16)
            for (int x = 0; x < frame.width; ++x) {
                // Calculate panning
                float pan; // Panning value [-1.0 (left) to 1.0 (right)]

                // Linear panning based on column position
                if ((frame.width > 1))
                    pan = 2.0f * (static_cast<float>(x) / static_cast<float>(frame.width - 1)) - 1.0f;
                else
                    pan = 0.0f;

                // Equal Power Panning
                float panL = 1.0f, panR = 1.0f; // Default für Mono
                if (useStereo) {
                    panL = std::sqrt(0.5f * (1.0f - pan)); // Left channel
                    panR = std::sqrt(0.5f * (1.0f + pan)); // Right channel
                }

                // Calculate amplitudes for each row based on pixel values
                for (int y = 0; y < frame.height; ++y) {
                    float c = frame.pixels[y * frame.width + x]; // Pixel value at (x,y)
                    if (params.scaleType == AudioParams::ScaleType::LINEAR || img.frames.size() > 1 /* Recognise GIF */) {
                        // old Sign-Pow Logic for LINEAR and GIFs
                        c -= 0.5f; // Center around 0
                        float sign = (c >= 0.0f) ? 1.0f : -1.0f;    // Sign of c
                        float ac = std::abs(c);                     // Absolute value of c
                        int li = std::min(static_cast<int>(ac * static_cast<float>(LUT_SIZE - 1)), LUT_SIZE - 1); // Lookup index
                        amp[y] = sign * gammaLUT[li]; // Apply gamma correction using LUT
                    } else {
                        // LOG / MEL / BARK für JPGs
                        c = std::max(c, 1e-4f);             // Avoid log(0)
                        amp[y] = std::pow(c, params.gamma); // Apply gamma correction
                    }

                    float startPhase = phaseInc[y] * static_cast<float>(x) * static_cast<float>(colSamples); // Start phase for this column
                    sine_values[y] = sinf(startPhase);      // Initial sine value
                    cosine_values[y] = cosf(startPhase);    // Initial cosine value
                }


                // Generate samples for this column
                const int colSamplesFrame = colSamples; // Samples per column for this frame (to avoid recalculating every time)
                for (int i = 0; i < colSamplesFrame; ++i) {
                    float sample = 0.0f;                        // Sample for this point in time
                    for (int y = 0; y < frame.height; ++y) {    // add up all frequencies
                        sample += amp[y] * sine_values[y];      // Add sine wave contribution

                        // Update sine and cosine using recursive formula
                        float s = sine_values[y];                           // Current sine value
                        float c = cosine_values[y];                         // Current cosine value
                        sine_values[y] = s * cosInc[y] + c * sinInc[y];     // Update sine
                        cosine_values[y] = c * cosInc[y] - s * sinInc[y];   // Update cosine
                    }
                    if (useWindow) sample *= hann[i]; // Apply windowing if needed

                    int idx = offset + x * colSamples + i; // Index in final audio
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

    // --- Remove DC ---

    // Remove DC-Offset
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

    // --- Normalise ---
    if (useStereo) { // Stereo
        if (norm == NormType::RMS) {    // RMS
            double sumSq = 0.0;         // Sum of squares
            
            // Add up squares of both channels to calculate RMS
            for (size_t i = 0; i < finalAudioL.size(); ++i) {
                sumSq += finalAudioL[i] * finalAudioL[i];
                sumSq += finalAudioR[i] * finalAudioR[i];
            }
           
            // RMS over both channels
            auto rms = static_cast<float>(std::sqrt(sumSq / (2.0 * static_cast<float>(finalAudioL.size()))));
            float gain = 0.1f / std::max(rms, 1e-6f);           // Calculate gain
            for (size_t i = 0; i < finalAudioL.size(); ++i) {   // Normalize both channels
                finalAudioL[i] *= gain;
                finalAudioR[i] *= gain;
            }
        } else {
            // PEAK
            float maxAmp = 0.0f; // Maximum amplitude
            // Find maximum absolute amplitude across both channels
            for (size_t i = 0; i < finalAudioL.size(); ++i) // Check both channels
                maxAmp = std::max({maxAmp, std::abs(finalAudioL[i]), std::abs(finalAudioR[i])});
            if (maxAmp > 0.0f) { // Normalise if maxAmp > 0
                float inv = 1.0f / maxAmp;                          // Inverse of max amplitude
                for (size_t i = 0; i < finalAudioL.size(); ++i) {   // Normalize both channels
                    finalAudioL[i] *= inv;
                    finalAudioR[i] *= inv;
                }
            }
        }
    } else {
        if (norm == NormType::RMS) {
            double sumSq = 0.0;
            for (float s: finalAudio) sumSq += s * s;
            auto rms = static_cast<float>(std::sqrt(sumSq / static_cast<float>(finalAudio.size()))); // Calculate RMS
            float gain = 0.1f / std::max(rms, 1e-6f); // Calculate gain
            for (float &s: finalAudio) s *= gain; // Normalize
        } else {
            // PEAK
            // Find maximum absolute amplitude
            float maxAmp = *ranges::max_element(finalAudio, [](const float a, const float b) -> bool {
                return std::abs(a) < std::abs(b);
            });
            if (maxAmp > 0.0f) { // Normalise if maxAmp > 0
                float inv = 1.0f / maxAmp;
                for (auto &s: finalAudio) s *= inv;
            }
        }
    }


    clog << "Audio post-processing completed." << endl;
    clog << "Saving audio to file..." << endl;

    // --- Save WAV ---
    // Save WAV
    if (!save_wav_file(useStereo, wavPath.string(), finalAudio, finalAudioL, finalAudioR, params.samplerate)) {
        cerr << "Could not save WAV file to " << wavPath.string() << endl;
        return 1;
    }
    clog << "WAV file saved successfully." << endl;

    // --- If needed convert to other format with Ffmpeg ---
    if (outputFormat != "wav") {
        if (!convert_file(outputFormat, wavPath.string(), keepWav, finalOutputPath)) { // Convert file
            cerr << "Could not convert file " << wavPath.string() << endl;
            cerr << "Wav file is kept." << endl;
            return 1;
        }
        clog << "File conversion completed Successfully." << endl;
    }

    clog << "Good bye! :3" << endl;
    return 0;
}

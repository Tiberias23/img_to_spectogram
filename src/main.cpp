#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <fstream>
#include <numeric>

#include <AudioFile.h>
#include <cxxopts.hpp>
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

using namespace std;

// --- Structs ---
/**
 * @brief Represents a single image frame with pixel data.
 */
struct ImageFrame {
    int width{};
    int height{};
    std::vector<float> pixels; // normalisierte Pixelwerte
    std::vector<int> delays; // nur für GIF
};

/**
 * @brief Simple image structure to hold frames.
 */
struct Image {
    int channels = 1; // Graustufen
    std::vector<ImageFrame> frames;

    /**
     * @brief load an image or GIF from the given path
     * @param path the path to load the image or GIF from
     * @return true if successful, false otherwise
     */
    bool load(const std::string &path) {
        if (path.find(".gif") != std::string::npos) {
            return loadGIF(path);
        }
        return loadImage(path);
    }

private:
    /**
     * @brief Loads a file into a vector of unsigned char.
     * @param path the path to load the file from
     * @return a vector with the frames
     * @author Lupo
     */
    static std::vector<unsigned char> load_file(const std::string &path) {
        std::ifstream file(path, std::ios::binary);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + path);
        }
        return std::vector<unsigned char>{
            std::istreambuf_iterator<char>(file),
            std::istreambuf_iterator<char>()
        };
    }

    /**
     * @brief Loads a GIF image from the specified path.
     * @param path The file path of the GIF image.
     * @return True if the GIF was loaded successfully, false otherwise.
     */
    bool loadGIF(const std::string &path) {
        const auto data = load_file(path);

        int *delays = nullptr;
        int frames_count = 0;
        int width = 0, height = 0;

        unsigned char *gif = stbi_load_gif_from_memory(
            data.data(), static_cast<int>(data.size()),
            &delays, &width, &height, &frames_count, nullptr, channels
        );

        if (!gif || frames_count <= 0) {
            std::cerr << "Failed to load gif '" << path << "': " << stbi_failure_reason() << std::endl;
            return false;
        }
        clog << "loaded gif '" << path << "' (" << width << "x" << height << "), " << frames_count << " frames" << endl;

        for (int f = 0; f < frames_count; ++f) {
            ImageFrame frame;
            frame.width = width;
            frame.height = height;
            frame.pixels.resize(width * height);
            frame.delays.resize(1); // ein Delay pro Frame

            for (int i = 0; i < width * height; ++i)
                frame.pixels[i] = static_cast<float>(gif[f * width * height + i]) / 255.0f;

            frame.delays[0] = (delays != nullptr) ? delays[f] : 10; // default 100ms
            frames.push_back(frame);
        }

        stbi_image_free(gif);
        if (delays) free(delays);

        return true;
    }

    /**
     * @brief Loads an image from the specified path.
     * @param path The file path of the image.
     * @return True if the image was loaded successfully, false otherwise.
     */
    bool loadImage(const std::string &path) {
        int w, h, c;
        unsigned char *data = stbi_load(path.c_str(), &w, &h, &c, channels);
        if (!data) return false;
        clog << "loaded image '" << path << "' (" << w << "x" << h << ")" << endl;
        ImageFrame frame;
        frame.width = w;
        frame.height = h;
        frame.pixels.resize(w * h);
        for (int i = 0; i < w * h; ++i)
            frame.pixels[i] = static_cast<float>(data[i]) / 255.0f;

        frames.push_back(frame);
        stbi_image_free(data);
        return true;
    }
};

/**
 * @brief Audio parameters for sound generation.
 */
struct AudioParams {
    int samplerate = 44100;
    float minFreq = 400.0f;
    float maxFreq = 14000.0f;
    float durationPerColumn = 0.01f;
};

// --- Functions ---

/**
 * @brief Normalizes the audio buffer to the range [-1.0, 1.0].
 * @param audio The audio buffer to normalize.
 */
inline void normalizeAudio(std::vector<float> &audio) {
    float maxAmp = 0.0f;
    for (const float v: audio) maxAmp = std::max(maxAmp, std::abs(v));
    if (maxAmp > 0.0f) {
        for (float &v: audio) v /= maxAmp;
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
            ("duration-per-column", "Duration per image column in seconds",
             cxxopts::value<float>()->default_value("0.01"));

    const cxxopts::ParseResult parsed = options.parse(argc, argv);
    if (parsed.count("help")) {
        std::cout << options.help() << std::endl;
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
        for (int y = 0; y < frame.height; ++y)
            frequencies[y] = params.minFreq * pow(params.maxFreq / params.minFreq,
                                                  static_cast<float>(frame.height - 1 - y) / static_cast<float>(
                                                      frame.height - 1));

        // Jede Spalte bekommt mindestens 256 Samples für scharfes Spektrogramm
        int colSamples = max(256, samplesPerColumn);

        // Hann-Window vorberechnen (reduziert Spectral Leakage)
        vector<float> hann(colSamples);
        for (int i = 0; i < colSamples; ++i) {
            hann[i] = 0.5f * static_cast<float>(1.0f - cos(2.0f * M_PI * i / (colSamples - 1)));
        }

        for (int x = 0; x < frame.width; ++x) {
            for (int i = 0; i < colSamples; ++i) {
                float globalT = static_cast<float>(offset + x * colSamples + i) / static_cast<float>(params.samplerate);
                float sample = 0.0f;
                for (int y = 0; y < frame.height; ++y) {
                    float centered = frame.pixels[y * frame.width + x] - 0.5f;
                    sample += centered * static_cast<float>(sin(2.0f * M_PI * frequencies[y] * globalT));
                }
                sample *= hann[i];
                if (int indexInAudioFile = offset + x * colSamples + i; indexInAudioFile < finalAudio.size())
                    finalAudio[indexInAudioFile] += sample;
            }
        }
    }

    float mean = std::accumulate(finalAudio.begin(), finalAudio.end(), 0.0f)
                 / static_cast<float>(finalAudio.size());

    for (auto &s: finalAudio)
        s -= mean;

    // Normalisieren
    normalizeAudio(finalAudio);

    // --- WAV speichern ---
    AudioFile<float> wav;
    wav.setNumChannels(1);
    wav.setSampleRate(params.samplerate);
    wav.setNumSamplesPerChannel(static_cast<int>(finalAudio.size()));
    for (int i = 0; i < finalAudio.size(); ++i)
        wav.samples[0][i] = finalAudio[i];

    wav.save(outputSoundPath);
    cout << "Audio saved to " << outputSoundPath << endl;
    return 0;
}

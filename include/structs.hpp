//
// Created by lupo on 24.01.26.
//

#ifndef IMG_TO_SPECTROGRAM_STRUCTS_H
#define IMG_TO_SPECTROGRAM_STRUCTS_H

#include <vector>
#include <iostream>
#include <string>
#include <fstream>

/**
 * @brief Normalisierungs typen für Audio
 * @author Lupo
 */
enum class NormType { PEAK, RMS };

/**
 * @brief Represents a single image frame with pixel data.
 * @author Lupo
*/
struct ImageFrame {
    int width{};
    int height{};
    std::vector<float> pixels; // normalisierte Pixelwerte
    std::vector<int> delays; // nur für GIF
};

/**
 * @brief Simple image structure to hold frames.
 * @author Lupo
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
        std::clog << "loaded gif '" << path << "' (" << width << "x" << height << "), " << frames_count << " frames" << std::endl;

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
        std::clog << "loaded image '" << path << "' (" << w << "x" << h << ")" << std::endl;
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
 * @author Lupo
 */
struct AudioParams {
    enum class ScaleType {
        LINEAR,
        LOGARITHMIC,
        MEL,
        BARK
    };

    int samplerate = 44100;
    float minFreq = 400.0f;
    float maxFreq = 14000.0f;
    float durationPerColumn = 0.01f;
    float gamma{1.0f};
    ScaleType scaleType = ScaleType::LOGARITHMIC;
};

/**
 * @brief Stereo Normalization modes
 * @author Lupo
 */
enum class StereoNorm {
    LINKED, // beide Kanäle gemeinsam (empfohlen)
    INDEPENDENT // L und R getrennt
};

#endif //IMG_TO_SPECTROGRAM_STRUCTS_H
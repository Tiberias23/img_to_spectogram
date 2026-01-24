#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>

#include <AudioFile.h>
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

using namespace std;

// --- Structs ---
struct ImageFrame {
    int width;
    int height;
    std::vector<float> pixels; // normalisierte Pixelwerte
};

struct Image {
    int channels = 1; // Graustufen
    std::vector<ImageFrame> frames;

    bool loadSingleImage(const std::string& path) {
        int w, h, c;
        unsigned char* data = stbi_load(path.c_str(), &w, &h, &c, channels);
        if (!data) return false;

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

struct AudioParams {
    int samplerate = 44100;
    float minFreq = 100.0f;
    float maxFreq = 10000.0f;
    float durationPerColumn = 0.01f;
};

struct FrequencyTable {
    std::vector<float> frequencies;
    std::vector<float> sin_table;

    void generate(int height, int samplesPerColumn, float minFreq, float maxFreq, int samplerate) {
        frequencies.resize(height);
        for (int y = 0; y < height; ++y) {
            frequencies[y] = minFreq * std::pow(maxFreq / minFreq, static_cast<float>(height - 1 - y) / (height - 1));
        }

        sin_table.resize(height * samplesPerColumn);
        for (int y = 0; y < height; ++y) {
            float freq = frequencies[y];
            for (int i = 0; i < samplesPerColumn; ++i) {
                float time = static_cast<float>(i) / samplerate;
                sin_table[y * samplesPerColumn + i] = std::sin(2.0f * M_PI * freq * time);
            }
        }
    }
};

// --- Main ---
int main() {
    string inputImagePath = "/home/lupo/Bilder/Boy_Bi_kissermeme/Silly_Cat_Character.jpg";
    string outputSoundPath = "./bild.wav";

    // --- Image laden ---
    Image img;
    if (!img.loadSingleImage(inputImagePath)) {
        cout << "Failed to load image." << endl;
        return 1;
    }
    ImageFrame& frame = img.frames[0]; // erstes Frame nutzen
    cout << "Image loaded: " << frame.width << "x" << frame.height << endl;

    // --- Audio-Parameter ---
    AudioParams params;
    int samplesPerColumn = static_cast<int>(params.samplerate * params.durationPerColumn);
    int totalSamples = samplesPerColumn * frame.width;
    vector<float> audio(totalSamples, 0.0f);

    // --- Frequenzen / Sinus-Tabelle ---
    FrequencyTable table;
    table.generate(frame.height, samplesPerColumn, params.minFreq, params.maxFreq, params.samplerate);

    // --- Audio erzeugen ---
    for (int x = 0; x < frame.width; ++x) {
        for (int i = 0; i < samplesPerColumn; ++i) {
            int sampleIndex = x * samplesPerColumn + i;
            audio[sampleIndex] = 0.0f;

            for (int y = 0; y < frame.height; ++y) {
                float amp = frame.pixels[y * frame.width + x];
                audio[sampleIndex] += amp * table.sin_table[y * samplesPerColumn + i];
            }
        }
    }

    // --- Normalisieren ---
    float maxAmp = 0.0f;
    for (float v : audio) maxAmp = std::max(maxAmp, std::abs(v));
    if (maxAmp > 0.0f) {
        for (float& v : audio) v /= maxAmp;
    }

    // --- WAV speichern ---
    AudioFile<float> wav;
    wav.setNumChannels(1);
    wav.setSampleRate(params.samplerate);
    wav.setNumSamplesPerChannel(totalSamples);
    for (int i = 0; i < totalSamples; ++i)
        wav.samples[0][i] = audio[i];

    wav.save(outputSoundPath);
    cout << "Audio saved to " << outputSoundPath << endl;

    return 0;
}

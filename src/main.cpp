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

// Für meine structs
#include <structs.hpp>
using namespace std;

// --- Funktionen ---

/**
 * @brief macht einen Systemcall um zu prüfen ob ffmpeg im PATH ist
 * @details macht einen Systemcall "ffmpeg -version" und prüft den Rückgabewert.
  *         stdout und stderr werden verworfen.
  *         Wenn der Rückgabewert 0 ist, wurde ffmpeg gefunden.
 * @details Dies ist notwendig, um Audio in andere Formate als WAV zu konvertieren.
 * @details ffmpeg muss installiert und im PATH sein.
  *          Andernfalls wird nur WAV ausgegeben.
  *@details  Siehe https://ffmpeg.org/download.html für Installationsanweisungen.
  *          Unter Windows kann ffmpeg von https://ffmpeg.org/download.html#build-windows heruntergeladen werden.
  *          Der entpackte Ordner muss dann zum PATH hinzugefügt werden.
  *          Alternativ kann ffmpeg auch in das gleiche Verzeichnis wie die ausführbare Datei kopiert werden.
 * @return true, wenn ffmpeg im PATH gefunden wurde
 * @author Lupo
 */
[[nodiscard]] inline bool ffmpegExists() {
#ifdef _WIN32
    return std::system("ffmpeg -version > NUL 2>&1") == 0;
#else
    return std::system("ffmpeg -version > /dev/null 2>&1") == 0;
#endif
}

/**
 * @brief Konvertiert eine WAV Datei mit ffmpeg in ein anderes Format
 * @details Konvertiert die Datei wavPath in targetPath mit ffmpeg.
 *            stdout und stderr werden in eine temporäre Log-Datei umgeleitet.
 *            Im Fehlerfall wird der Inhalt der Log-Datei ausgegeben.
 *            Die Log-Datei wird im Erfolgsfall gelöscht.
 *            Dies geschieht, um das Konsolenfenster sauber zu halten. :3
 * @param wavPath the input wav file path
 * @param targetPath the output file path
 * @return true if conversion was successful false otherwise
 * @author Lupo
 */
bool convertWithFFmpeg(const std::string &wavPath, const std::string &targetPath) {
    // Temporäre Log-Datei
    std::string logFile = "ffmpeg.log";

    // ffmpeg-Aufruf: stdout und stderr in logFile umleiten
    std::string cmd = "ffmpeg -y -i \"" + wavPath + "\" \"" + targetPath + "\" > \"" + logFile + "\" 2>&1";
    const int ret = std::system(cmd.c_str());

    if (ret != 0) {
        // Im Fehlerfall Log-Datei lesen und ausgeben
        const std::ifstream log(logFile);
        if (log.is_open()) {
            std::stringstream buffer;
            buffer << log.rdbuf();
            std::cerr << "ffmpeg failed:\n" << buffer.str() << std::endl;
        } else { // Für den Fall, dass die log datei nicht erstellt werden oder gelesen werden konnte, dass sollte aber eigentlich nie passieren.
            std::cerr << "ffmpeg failed and log file could not be opened." << std::endl;
        }
        return false;
    }

    // Erfolg → Log löschen
    std::remove(logFile.c_str());
    return true;
}

/**
 * @brief Wandelt Hz in Mel-Frequenz um
 * @param f die Frequenz in Hz
 * @return die Mel-Frequenz
 * @author Lupo
 */
inline float hzToMel(const float f) {
    return 2595.0f * std::log10(1.0f + f / 700.0f);
}

/**
 * @brief Wandelt Mel-Frequenz in Hz um
 * @param m die Mel-Frequenz
 * @return die Frequenz in Hz
 * @author Lupo
 */
inline float melToHz(const float m) {
    return 700.0f * (std::pow(10.0f, m / 2595.0f) - 1.0f);
}

/**
 * @brief Wandelt Hz in Bark-Frequenz um
 * @param f die Frequenz in Hz
 * @return die Bark-Frequenz
 * @author Lupo
 */
inline float hzToBark(const float f) {
    return 13.0f * std::atan(0.00076f * f)
         + 3.5f * std::atan(std::pow(f / 7500.0f, 2.0f));
}

/**
 * @brief Wandelt Bark-Frequenz in Hz um
 * @param z die Bark-Frequenz
 * @return die Frequenz in Hz
 * @author Lupo
 */
inline float barkToHz(const float z) {
    // numerische Inversion (Newton-Raphson wäre overkill)
    // → wir benutzen eine gute Näherung
    return 600.0f * std::sinh(z / 6.0f);
}


/**
 * @brief Berechnet die Frequenz für eine gegebene y-Position im Bild
 * @param y die y-Position (0 = unten, height-1 = oben)
 * @param height die Höhe des Bildes
 * @param params die Audio-Parameter
 * @return die Frequenz in Hz
 * @author Lupo
 */
float freqForY(const int y, const int height, const AudioParams& params) {
    if (height <= 1) // Vermeidung Division durch Null
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
            ("norm", "Normalization: peak|rms", cxxopts::value<std::string>()->default_value("peak"));

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
    params.gamma = parsed["gamma"].as<float>();

    auto toLower = [](std::string s) {
        ranges::transform(s, s.begin(), ::tolower);
        return s;
    };

    std::string scaleStr = toLower(parsed["scale"].as<std::string>());
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

    std::string normStr = toLower(parsed["norm"].as<std::string>());
    enum class NormType { PEAK, RMS } norm;

    if (normStr == "peak") {
        norm = NormType::PEAK;
    } else if (normStr == "rms") {
        norm = NormType::RMS;
    } else {
        throw std::runtime_error("Invalid normalization type: " + normStr + " (use peak|rms)");
    }


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
    clog << "Generating audio: " << totalSamples << " samples at " << params.samplerate << " Hz" << endl;
    clog << "This may take a while depending on image/gif size and number of frames..." << endl;
    vector<float> finalAudio(totalSamples, 0.0f);
    for (size_t f = 0; f < img.frames.size(); ++f) {
        auto &frame = img.frames[f];
        int offset = frameOffsets[f];
        vector<float> frequencies(frame.height);
        vector<float> phaseInc(frame.height);
        for (int y = 0; y < frame.height; ++y) {
            frequencies[y] = freqForY(y, frame.height, params);             // Frequenz für diese Zeile
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
        #pragma omp parallel for schedule(static) default(none) \
            shared(finalAudio, frame, frequencies, hann, offset, colSamples,params, phaseInc)
        for (int x = 0; x < frame.width; ++x) {
            vector<float> phase(frame.height, 0.0f); // Phase pro Frequenz
            vector<float> amp(frame.height); // Amplitude pro Frequenz
            // Amplitude vorberechnen
            for (int y = 0; y < frame.height; ++y) {
                float c = frame.pixels[y * frame.width + x] - 0.5f;
                float sign = (c >= 0.0f) ? 1.0f : -1.0f;
                amp[y] = sign * std::pow(std::abs(c), params.gamma);
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
    clog << "Audio generation completed." << endl;
    clog << "Post-processing audio..." << endl;
    // --- Normalisieren ---

    float mean = std::accumulate(finalAudio.begin(), finalAudio.end(), 0.0f)
                 / static_cast<float>(finalAudio.size());

    for (auto &s: finalAudio)
        s -= mean;

    // Normalisieren
    if (norm == NormType::RMS) {
        float rms = 0.0f;
        for (float s : finalAudio) rms += s * s;
        rms = std::sqrt(rms / static_cast<float>(finalAudio.size()));

        float target = 0.1f;
        float gain = target / std::max(rms, 1e-6f);
        for (float& s : finalAudio) s *= gain;
    } else if (norm == NormType::PEAK) { // Peak-Normalisierung
        float maxAmp = 0.0f;
        for (const float v: finalAudio) maxAmp = std::max(maxAmp, std::abs(v));
        if (maxAmp > 0.0f) {
            for (float &v: finalAudio) v /= maxAmp;
        }
    }


    clog << "Audio post-processing completed." << endl;
    clog << "Saving audio to file..." << endl;

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
        clog << "Converting WAV to " << outputFormat << " using ffmpeg..." << endl;
        if (!ffmpegExists()) {
            cerr << "ffmpeg is not installed or not found in PATH. Cannot convert to " << outputFormat << endl;
            cerr << "Will keep the WAV file at " << outputSoundPath << endl;
            return 1;
        }
        if (!convertWithFFmpeg(outputSoundPath, finalOutputPath.string())) {
            cerr << "Failed to convert WAV to " << outputFormat << endl;
            cerr << "Make sure ffmpeg supports the format." << endl;
            cerr << "WAV file is kept at " << outputSoundPath << endl;
            cerr << "See ffmpeg.log for details." << endl;
            return 1;
        }
        // Original WAV löschen
        if (!keepWav) {
            std::remove(outputSoundPath.c_str());
            cout << "Intermediate WAV file deleted." << endl;
        }
        cout << "Converted audio saved to " << finalOutputPath.string() << endl;
    }
    clog << "Good bye! :3" << endl;
    return 0;
}

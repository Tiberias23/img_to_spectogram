//
// Created by lupo on 25.01.26.
//

#ifndef IMG_TO_SPECTROGRAM_UNGROUPED_UTILITY_FUNKTIONS_HPP
#define IMG_TO_SPECTROGRAM_UNGROUPED_UTILITY_FUNKTIONS_HPP
/**
 * @brief makes a Syscall to check if ffmpeg is installed
 * @details uses "ffmpeg -version" command to check if ffmpeg is installed
  *          and available in PATH.
  *          Redirects output to null device to keep console clean.
  *          Works on Windows and Unix-like systems.
 * @return true if ffmpeg is found in PATH false otherwise
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
 * @brief Convert a wav file to another format using ffmpeg
 * @details Converts a wav file to another format using ffmpeg.
 *            stdout and stderr are diverted to a log file.
 *            In case of failure the log file is read and printed to the console.
 *            The Logfile is removed if everything went well.
 *            this happens to have a clean console. :3
 * @param wavPath the input wav file path
 * @param targetPath the output file path
 * @return true if conversion was successful false otherwise
 * @author Lupo
 */
inline bool convertWithFFmpeg(const std::string &wavPath, const std::string &targetPath) {
    // Temporary log file for ffmpeg output
    std::string logFile = "ffmpeg.log";

    // ffmpeg command with redirected output to a log file
    std::string cmd = "ffmpeg -y -i \"" + wavPath + "\" \"" + targetPath + "\" > \"" + logFile + "\" 2>&1";
    const int ret = std::system(cmd.c_str());

    if (ret != 0) {
        // if ffmpeg failed → read and print log file
        const std::ifstream log(logFile);
        if (log.is_open()) {
            std::stringstream buffer;
            buffer << log.rdbuf();
            std::cerr << "ffmpeg failed:\n" << buffer.str() << std::endl;
        } else { // If ffmpeg failed and the log file could not be opened should never happen
            std::cerr << "ffmpeg failed and log file could not be opened." << std::endl;
        }
        return false;
    }

    // Success → Remove log file
    std::remove(logFile.c_str());
    return true;
}

/**
 * @brief Convert Hz to Mel-Frequenz
 * @param f the Frequency in Hz
 * @return the Mel-Frequency
 * @author Lupo
 */
inline float hzToMel(const float f) {
    return 2595.0f * std::log10(1.0f + f / 700.0f);
}

/**
 * @brief convert Mel-Frequenz in Hz
 * @param m the Mel-Frequency
 * @return the Frequency in Hz
 * @author Lupo
 */
inline float melToHz(const float m) {
    return 700.0f * (std::pow(10.0f, m / 2595.0f) - 1.0f);
}

/**
 * @brief Convert Hz to Bark-Frequenz
 * @param f the Frequency in Hz
 * @return the Bark-Frequency
 * @author Lupo
 */
inline float hzToBark(const float f) {
    return 13.0f * std::atan(0.00076f * f)
         + 3.5f * std::atan(std::pow(f / 7500.0f, 2.0f));
}

/**
 * @brief Convert Bark-Frequenz to Hz
 * @param z the Bark-Frequency
 * @return the Frequency in Hz
 * @author Lupo
 */
inline float barkToHz(const float z) {
    // numerische Inversion (Newton-Raphson wäre overkill)
    // → we use approximation formula
    return 600.0f * std::sinh(z / 6.0f);
}

/**
 * @brief Converts a WAV file to another format using ffmpeg
 * @param outputFormat the desired output format e.g. "mp3", "flac", "wav"
 * @param outputSoundPath the path to the intermediate wav file
 * @param keepWav whether to keep the intermediate wav file
 * @param finalOutputPath the final output file path (has to be set if outputFormat is not "wav")
 * @return true if conversion was successful false otherwise
 * @author Lupo
 */
[[nodiscard]] inline bool convert_file(const std::string &outputFormat, const std::string &outputSoundPath, const bool keepWav,
                                              const std::filesystem::path &finalOutputPath = "./" /*Dummy Default*/) {
    std::clog << "Converting WAV to " << outputFormat << " using ffmpeg..." << std::endl;
    if (!ffmpegExists()) {
        std::cerr << "ffmpeg is not installed or not found in PATH. Cannot convert to " << outputFormat << std::endl;
        std::cerr << "Will keep the WAV file at " << outputSoundPath << std::endl;
        return false;
    }
    if (!convertWithFFmpeg(outputSoundPath, finalOutputPath.string())) {
        std::cerr << "Failed to convert WAV to " << outputFormat << std::endl;
        std::cerr << "Make sure ffmpeg supports the format." << std::endl;
        std::cerr << "WAV file is kept at " << outputSoundPath << std::endl;
        std::cerr << "See ffmpeg.log for details." << std::endl;
        return false;
    }
    // Remove the intermediate WAV file if not needed
    if (!keepWav) {
        std::remove(outputSoundPath.c_str());
        std::cout << "Intermediate WAV file deleted." << std::endl;
    }
    std::cout << "Converted audio saved to " << finalOutputPath.string() << std::endl;
    return true;
}

/**
 * @brief Calculates the frequency for a given y position in the image
  *        according to the selected scale type in AudioParams.
 * @param y the y position in the image
 * @param height the height of the image
 * @param params the audio parameters
 * @return the frequency in Hz for the given y position
 * @author Lupo
 */
inline float freqForY(const int y, const int height, const AudioParams &params) {
    if (height <= 1) // Avoid division with zero
        return params.minFreq;

    const float t = static_cast<float>(height - 1 - y) / static_cast<float>(height - 1); // Normalized position (0 at bottom, 1 at top)

    switch (params.scaleType) {
        case AudioParams::ScaleType::LINEAR:
            return params.minFreq + t * (params.maxFreq - params.minFreq);

        case AudioParams::ScaleType::LOGARITHMIC:
            return params.minFreq * std::pow(params.maxFreq / params.minFreq, t);

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
    return params.minFreq; // unreachable, but the Compiler is happy
}

#endif //IMG_TO_SPECTROGRAM_UNGROUPED_UTILITY_FUNKTIONS_HPP
//
// Created by lupo on 25.01.26.
//

#ifndef IMG_TO_SPECTROGRAM_UNGROUPED_UTILITY_FUNKTIONS_HPP
#define IMG_TO_SPECTROGRAM_UNGROUPED_UTILITY_FUNKTIONS_HPP
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
inline bool convertWithFFmpeg(const std::string &wavPath, const std::string &targetPath) {
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
 * @brief Konvertiert die WAV-Datei in ein anderes Format mit ffmpeg
 * @param outputFormat das gewünschte Ausgabeformat (z.B. "mp3", "flac")
 * @param outputSoundPath der Pfad zur WAV-Datei
 * @param keepWav ob die WAV-Datei behalten werden soll
 * @param finalOutputPath der Pfad zur finalen Ausgabedatei (standardmäßig "./" als Dummy-Parameter, der muss übergeben werden!!!)
 * @return true bei Erfolg, false bei Fehler
 * @author Lupo
 */
[[nodiscard]] inline bool convert_file(const std::string &outputFormat, const std::string &outputSoundPath, const bool keepWav,
                                              const std::filesystem::path &finalOutputPath = "./" /*Dummy Parameter*/) {
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
    // Original WAV löschen
    if (!keepWav) {
        std::remove(outputSoundPath.c_str());
        std::cout << "Intermediate WAV file deleted." << std::endl;
    }
    std::cout << "Converted audio saved to " << finalOutputPath.string() << std::endl;
    return true;
}

/**
 * @brief Berechnet die Frequenz für eine gegebene y-Position im Bild
 * @param y die y-Position (0 = unten, height-1 = oben)
 * @param height die Höhe des Bildes
 * @param params die Audio-Parameter
 * @return die Frequenz in Hz
 * @author Lupo
 */
inline float freqForY(const int y, const int height, const AudioParams &params) {
    if (height <= 1) // Vermeidung Division durch null
        return params.minFreq;

    const float t = static_cast<float>(height - 1 - y) /
                    static_cast<float>(height - 1);

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
    return params.minFreq; // unreachable, aber Compiler happy
}

#endif //IMG_TO_SPECTROGRAM_UNGROUPED_UTILITY_FUNKTIONS_HPP
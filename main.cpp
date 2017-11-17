#include <QDateTime>
#include <QFile>
#include <QStringList>
#include <QTextStream>

#include <cstdint>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>
#include <cmath>

#include "immintrin.h"

// declare MSSIM function, which is implemented at the end
double mssim_11x11(const unsigned char *image1, const unsigned char *image2, const int imageWidth, const int imageHeight);

// all logic is in the main function including parallelism through lambda subroutines
int main(int argc, char *argv[])
{
    const uint64_t appStartTime = QDateTime::currentMSecsSinceEpoch();

    QString appName = "video-quality-metrics";
    QString msgPrefix = "[" + appName + "] ";

    // ---------------------------------------------------------------------------------------------- //
    // -- parse command line arguments -------------------------------------------------------------- //
    // ---------------------------------------------------------------------------------------------- //

    // check if we have any arguments at all
    if (argc == 1) {
        QTextStream(stderr, QIODevice::WriteOnly) << msgPrefix << "error: no arguments" << '\n';
        return -1;
    }

    // essential arguemnts (no need for separate flags as strings can be empty)
    QString videoFilePath1 = "";
    QString videoFilePath2 = "";

    // optional arguements
    bool noPSNR = false;
    bool noMSSIM = false;
    bool showProgress = true;
    bool showTime = true;
    bool printPerFrame = false;
    uint32_t numThreads = 1;

    // for each arguemnt
    for (int i = 1;  i < argc;  i++)
    {
        QString currentArgument (argv[i]);

        if (currentArgument == "-h" || currentArgument == "--help") {
            QTextStream(stdout, QIODevice::WriteOnly) << "syntax: " << appName << " <video_file_1> <video_file_2> [-nopsnr|-nomssim] [-noprogress] [-notime] [-print-per-frame] [-threads <int>]" << '\n';
            return 0;
        } else if (currentArgument == "-nopsnr") {
            noPSNR = true;
        } else if (currentArgument == "-nomssim") {
            noMSSIM = true;
        } else if (currentArgument == "-noprogress") {
            showProgress = false;
        } else if (currentArgument == "-notime") {
            showTime = false;
        } else if (currentArgument == "-print-per-frame") {
            printPerFrame = true;
        } else if (currentArgument == "-threads") {
            i++;
            if (i == argc) {
                QTextStream(stderr, QIODevice::WriteOnly) << msgPrefix << "error: missing value after \"-threads\"" << '\n';
                return -1;
            }
            bool ok;
            numThreads = QString(argv[i]).toUInt(&ok);
            if (!ok || numThreads == 0) {
                QTextStream(stderr, QIODevice::WriteOnly) << msgPrefix << "error: incorrect value after \"-threads\"" << '\n';
                return -1;
            }
        } else {
            // first and second "unidentified" arguments are considered as names of input video files;
            // anything beyond that will trigger error
            if (videoFilePath1.isEmpty()) {
                videoFilePath1 = currentArgument;
            } else if (videoFilePath2.isEmpty()) {
                videoFilePath2 = currentArgument;
            } else {
                QTextStream(stderr, QIODevice::WriteOnly) << msgPrefix << QString("error: unknown option \"%1\"").arg(currentArgument) << '\n';
                return -1;
            }
        }
    }

    // check if all arguments are in place
    {
        if (videoFilePath1.isEmpty()) {
            QTextStream(stderr, QIODevice::WriteOnly) << msgPrefix << "error: missing input video files" << '\n';
            return -1;
        }

        if (videoFilePath2.isEmpty()) {
            QTextStream(stderr, QIODevice::WriteOnly) << msgPrefix << "error: missing second video file" << '\n';
            return -1;
        }

        if (noPSNR && noMSSIM) {
            QTextStream(stderr, QIODevice::WriteOnly) << msgPrefix << "error: -nopsnr and -nossim can't be used simultaneously" << '\n';
            return -1;
        }
    }

    // ---------------------------------------------------------------------------------------------- //
    // -- try to open input video files ------------------------------------------------------------- //
    // ---------------------------------------------------------------------------------------------- //

    QFile videoFile1 (videoFilePath1);

    if (!videoFile1.open(QFile::ReadOnly)) {
        QTextStream(stderr, QIODevice::WriteOnly) << msgPrefix << QString("error: can't open file \"%1\"").arg(videoFilePath1) << '\n';
        return -1;
    }

    QFile videoFile2 (videoFilePath2);

    if (!videoFile2.open(QFile::ReadOnly)) {
        QTextStream(stderr, QIODevice::WriteOnly) << msgPrefix << QString("error: can't open file \"%1\"").arg(videoFilePath2) << '\n';
        return -1;
    }

    // ---------------------------------------------------------------------------------------------- //
    // -- define function, which should get frame dimensions and video length ----------------------- //
    // -- just by parsing file header --------------------------------------------------------------- //
    // ---------------------------------------------------------------------------------------------- //

    auto getVideoDimensions = [msgPrefix] (QFile &inputFile, uint32_t &width, uint32_t &height, uint32_t &length) -> bool
    {
        // read header byte by byte
        QString header;

        // find where first frame starts
        {
            char c;
            do {
                inputFile.read(&c, 1);
                header.append(QChar(c));
            } while (c != '\n');
        }

        // remove trailing '\n'
        header.chop(1);

        // split header by space
        QStringList headerChunks = header.split(QChar(' '));

        // chech header magick
        if (headerChunks.at(0) != "YUV4MPEG2") {
            QTextStream(stderr, QIODevice::WriteOnly) << msgPrefix << "error: wrong header magick in file " << inputFile.fileName() << '\n';
            return false;
        }

        // extract frame dimensions
        width = 0;
        height = 0;

        for (int i = 0;  i < headerChunks.size();  i++)
        {
            QString currentChunk = headerChunks.at(i);

            if (currentChunk.startsWith(QChar('W'))) {
                width = currentChunk.remove(0, 1).toUInt();
            } else if (currentChunk.startsWith(QChar('H'))) {
                height = currentChunk.remove(0, 1).toUInt();
            } else if (currentChunk.startsWith(QChar('F'))) {
                //
            } else if (currentChunk.startsWith(QChar('I'))) {
                //
            } else if (currentChunk.startsWith(QChar('A'))) {
                //
            } else if (currentChunk.startsWith(QChar('C'))) {
                if (currentChunk != "C420jpeg" && currentChunk != "C420mpeg2" && currentChunk != "C420") {
                    QTextStream(stderr, QIODevice::WriteOnly) << msgPrefix << "error: only 4:2:0 colour subsampling ('420jpeg', '420mpeg2' and '420') is supported by this program" << '\n';
                    return false;
                }
            }
        }

        if (width == 0 || height == 0) {
            QTextStream(stderr, QIODevice::WriteOnly) << msgPrefix << "error: wrong video resolution in video file " << inputFile.fileName() << '\n';
            return false;
        }

        // count expected number of frames
        const uint64_t contentSize = inputFile.size() - inputFile.pos();
        const uint64_t yuvFrameSize = width * height * 3 / 2;
        const uint64_t frameSizeWithHeader = 6 + yuvFrameSize;

        length = contentSize / frameSizeWithHeader;

        // a sanity check
        if (contentSize % frameSizeWithHeader != 0) {
            QTextStream(stderr, QIODevice::WriteOnly) << msgPrefix << "error: file size doesn't correspond to the whole number of frames in file " << inputFile.fileName() << '\n';
            return false;
        }

        return true;
    };

    // ---------------------------------------------------------------------------------------------- //
    // -- obtain and compare video dimensions; ------------------------------------------------------ //
    // -- progrm is terminates immediately if error occures when parsing the header ----------------- //
    // ---------------------------------------------------------------------------------------------- //

    uint32_t width1;
    uint32_t height1;
    uint32_t length1;
    if (getVideoDimensions(videoFile1, width1, height1, length1) == false) return -1;

    uint32_t width2;
    uint32_t height2;
    uint32_t length2;
    if (getVideoDimensions(videoFile2, width2, height2, length2) == false) return -1;

    if (width1 != width2 || height1 != height2) {
        QTextStream(stderr, QIODevice::WriteOnly) << msgPrefix << QString("error: different frame dimensions (%1x%2 and %3x%4)").arg(width1).arg(height1).arg(width2).arg(height2) << '\n';
        return -1;
    }

    if (length1 != length2) {
        QTextStream(stderr, QIODevice::WriteOnly) << msgPrefix << QString("error: different number of frames (%1 and %2)").arg(length1).arg(length2) << '\n';
        return -1;
    }

    // ---------------------------------------------------------------------------------------------- //
    // -- define function for reading a single frame ------------------------------------------------ //
    // ---------------------------------------------------------------------------------------------- //

    // make some constants explicit
    const uint32_t frameWidth = width1;
    const uint32_t frameHeight = height1;
    const uint32_t numFrames = length1;

    const uint32_t yComponentSize = frameWidth * frameHeight;
    const uint32_t uvComponentSize = yComponentSize / 2;

    auto getNextFrame = [msgPrefix, yComponentSize, uvComponentSize] (QFile *videoFile, uint8_t *yChannel) -> bool
    {
        // assuming that frames are read sequentially,
        // we need to check only the frame header
        const char frameHeader [] = "FRAME\n";
        char tempBuffer [6];
        videoFile->read(tempBuffer, 6);
        for (int i = 0;  i < 6;  i++) {
            if (frameHeader[i] != tempBuffer[i]) {
                QTextStream(stderr, QIODevice::WriteOnly) << msgPrefix << "error: missing frame header in " << videoFile->fileName() << '\n';
                return false;
            }
        }

        // read only Y channel and skip the UV components
        videoFile->read((char *) yChannel, yComponentSize);
        videoFile->seek(videoFile->pos() + uvComponentSize);
        return true;
    };

    // ---------------------------------------------------------------------------------------------- //
    // -- define function, which is called from each thread to obtain next pair of frames ----------- //
    // ---------------------------------------------------------------------------------------------- //

    // frame counter; numeration starts from 1
    uint32_t currentFrame = 1;

    // mutex guarding entrance to this function
    std::mutex frameRequestMutex;

    auto getFramePair = [&videoFile1, &videoFile2, numFrames, getNextFrame, &currentFrame, &frameRequestMutex] (uint8_t *yChannel1, uint8_t *yChannel2) -> uint32_t
    {
        // only one thread at a time enters this function
        frameRequestMutex.lock();

        // check if we've run out of frames
        if (currentFrame > numFrames) {
            // in such case return invalid frame number
            frameRequestMutex.unlock();
            return 0xffffffff;
        } else {
            // read frames from both videos
            const bool ok1 = getNextFrame(&videoFile1, yChannel1);
            const bool ok2 = getNextFrame(&videoFile2, yChannel2);

            // return error code if something went wrong
            if ((ok1 && ok2) == false) {
                frameRequestMutex.unlock();
                return 0xffffffff;
            }

            // return current frame number and increase counter
            const uint32_t frameNumber = currentFrame;
            currentFrame++;
            frameRequestMutex.unlock();
            return frameNumber;
        }
    };

    // ---------------------------------------------------------------------------------------------- //
    // -- prepare function, which performs deterministic accumulation of the results ---------------- //
    // -- and prints information to the console if necessary ---------------------------------------- //
    // ---------------------------------------------------------------------------------------------- //

    // variables for target values
    uint64_t total_sum_of_squared_errors = 0;
    double total_sum_of_frame_mssim = 0;

    // precision for printing (digits after comma)
    const int psnrPrecision = 5;
    const int mssimPrecision = 7;

    // prepare command line output
    QTextStream cout (stdout, QIODevice::WriteOnly);

    // print header in case of per frame output
    if (printPerFrame) {
        cout << "Frame\t" << (noPSNR ? "" : "Y-PSNR\t\t") << (noMSSIM ? "" :"Y-MSSIM") << '\n';
    }

    // last processed frame number (must be initialized with the previous number to start frame)
    volatile uint32_t lastProcessedFrame = 0;

    // mutex used for deterministic addition (and printing to the console)
    std::mutex accumulatorMutex;

    // wait condition triggered every time quality metrics are added to accumulators
    std::condition_variable frameProcessedCondition;

    auto submitResult = [noPSNR, noMSSIM, showProgress, printPerFrame, yComponentSize, numFrames,
                         &total_sum_of_squared_errors, &total_sum_of_frame_mssim, psnrPrecision, mssimPrecision,
                         &cout, &lastProcessedFrame, &accumulatorMutex, &frameProcessedCondition] (uint32_t frameNumber, uint64_t frameSumOfSqrErrors, double frameMSSIM)
    {
        // make only one thread enter this function
        // by locking function access mutex with special wrapper 'unique_lock'
        std::unique_lock<std::mutex> lock(accumulatorMutex);

        // label to facilitate logic of awoken threads (it can be replaced with infinite loop)
        TRY_ACCESS_ACCUMULATORS:

        // check if current frame number is exactly after last one, which metrics were added to accumulators
        if (frameNumber == (lastProcessedFrame + 1))
        {
            // add metrics to accumulators
            total_sum_of_squared_errors += frameSumOfSqrErrors;
            total_sum_of_frame_mssim += frameMSSIM;

            // print notifications if needed (in separate scope for clarity)
            {
                if (printPerFrame)
                {
                    cout << frameNumber << '\t';

                    if (!noPSNR) {
                        const double framePSNR = 10 * log10((uint64_t) 255 * 255 * yComponentSize / (double) frameSumOfSqrErrors);
                        cout << QString::number(framePSNR, 'f', psnrPrecision) << '\t';
                    }

                    if (!noMSSIM) {
                        cout << QString::number(frameMSSIM, 'f', mssimPrecision);
                    }

                    cout << '\n';
                    cout.flush();
                }
                else if (showProgress)
                {
                    cout << '\r' << frameNumber << '/' << numFrames;
                    cout.flush();
                }
            }

            // increase number of processed frames
            lastProcessedFrame++;

            // notify any other thread, which may be waiting
            frameProcessedCondition.notify_all();
        }
        else
        {
            // otherwise wait for other threads to add their results to accumulators
            frameProcessedCondition.wait(lock);

            // attempt write to the file (only if condition is met, spurious wakeups are avoided with "if" statement)
            goto TRY_ACCESS_ACCUMULATORS;
        }
    };

    // ---------------------------------------------------------------------------------------------- //
    // -- define function, which is executed by each thread and performs calculation of ------------- //
    // -- sum of squared errors and Y-MSSIM for a pair of frames ------------------------------------ //
    // ---------------------------------------------------------------------------------------------- //

    auto compareFrames = [noPSNR, noMSSIM, frameWidth, frameHeight, yComponentSize, getFramePair, submitResult] ()
    {
        // prepare Y buffers for both frames
        uint8_t *image1 = new uint8_t [yComponentSize];
        uint8_t *image2 = new uint8_t [yComponentSize];

        // infinite loop over pairs of frames
        while (true)
        {
            // get next pair of frames
            const uint32_t frameNumber = getFramePair(image1, image2);
            if (frameNumber == 0xffffffff) break;

            // calculate sum of squared errors for Y-PSNR
            // initialize sum regardless of the fact if it's needed
            uint64_t frameSumOfSqrErrors = 0;

            // do nothing if PSNR is not needed
            if (!noPSNR) {
                for (uint32_t i = 0;  i < yComponentSize;  i++) {
                    const int32_t difference = (int32_t) image1[i] - (int32_t) image2[i];
                    frameSumOfSqrErrors += difference * difference;
                }
            }

            // calculate MSSIM if necessary
            const double frameMSSIM = noMSSIM ? 0 : mssim_11x11(image1, image2, frameWidth, frameHeight);

            // save results
            submitResult(frameNumber, frameSumOfSqrErrors, frameMSSIM);
        }

        // cleaning
        delete [] image1;
        delete [] image2;
    };

    // ---------------------------------------------------------------------------------------------- //
    // -- actual parallel processing ---------------------------------------------------------------- //
    // ---------------------------------------------------------------------------------------------- //

    std::vector<std::thread> threads (numThreads);
    for (uint32_t i = 0;  i < numThreads;  i++) threads[i] = std::thread(compareFrames);
    for (uint32_t i = 0;  i < numThreads;  i++) threads[i].join();

    // ---------------------------------------------------------------------------------------------- //
    // -- report time ------------------------------------------------------------------------------- //
    // ---------------------------------------------------------------------------------------------- //

    // close input files
    videoFile1.close();
    videoFile2.close();

    // remove dynamic line and erase it with spaces
    if (showProgress && !printPerFrame) {
        cout << '\r';
        for (uint32_t i = 0;  i < 100;  i++) cout << ' ';
        cout << '\r';
    }

    // finish calculating quality metrics
    const double psnr = 10 * log10((uint64_t) 255 * 255 * yComponentSize * numFrames / (double) total_sum_of_squared_errors);
    const double mssim = total_sum_of_frame_mssim / numFrames;

    // prepare metrics for printing
    QString resultString = "average ";
    resultString        += noPSNR  ? "" : (QString::number(psnr, 'f', psnrPrecision) + ',');
    resultString        += noMSSIM ? "" : (QString::number(mssim, 'f', mssimPrecision) + ',');

    // calculate elapsed time if necessary
    if (showTime) {
        const double elapsedTime = (QDateTime::currentMSecsSinceEpoch() - appStartTime) / 1000.0;
        resultString += QString::number(elapsedTime, 'f', 3) + ',';
    }

    // remove last comma (regardless of how many values need to be printed)
    resultString.chop(1);

    // print results
    cout << resultString << '\n';

    return 0;
}



// ---------------------------------------------------------------------------------------------- //
// -- this functions calculates mean SSIM between two grayscale images -------------------------- //
// -- using 11x11 sliding window moving with 2 pixels step horizontally and vertically ---------- //
// ---------------------------------------------------------------------------------------------- //

double mssim_11x11(const unsigned char *image1, const unsigned char *image2, const int imageWidth, const int imageHeight)
{
    // vector data type for AVX instructions
    typedef __m256 fp32x8;

    // normalized Gauss weighting coefficients for 11x11 pixels window (sigma of the distribution is 1.5)
    // extra 7 zeroes at the end are needed to fit a whole number of fp32x8 vectors (128 / 8 = 16 vectors)
    const float gaussNormalizedWeights [128] __attribute__ ((aligned (sizeof(fp32x8)))) =
    {
        1.05756557075783e-06, 7.8144113306272e-06, 3.70224761237097e-05, 0.000112464352203375, 0.000219050647191676, 0.000273561152999412, 0.000219050647191676, 0.000112464352203375, 3.70224761237097e-05, 7.8144113306272e-06, 1.05756557075783e-06,
        7.8144113306272e-06, 5.77411237021237e-05, 0.000273561152999412, 0.000831005407560634, 0.00161857752060636, 0.0020213587060008, 0.00161857752060636, 0.000831005407560634, 0.000273561152999412, 5.77411237021237e-05, 7.8144113306272e-06,
        3.70224761237097e-05, 0.000273561152999412, 0.00129605556026987, 0.00393706916085999, 0.00766836362659358, 0.00957662724216503, 0.00766836362659358, 0.00393706916085999, 0.00129605556026987, 0.000273561152999412, 3.70224761237097e-05,
        0.000112464352203375, 0.000831005407560634, 0.00393706916085999, 0.0119597601002285, 0.0232944318700625, 0.0290912248949643, 0.0232944318700625, 0.0119597601002285, 0.00393706916085999, 0.000831005407560634, 0.000112464352203375,
        0.000219050647191676, 0.00161857752060636, 0.00766836362659358, 0.0232944318700625, 0.0453713579203497, 0.0566619690238992, 0.0453713579203497, 0.0232944318700625, 0.00766836362659358, 0.00161857752060636, 0.000219050647191676,
        0.000273561152999412, 0.0020213587060008, 0.00957662724216503, 0.0290912248949643, 0.0566619690238992, 0.0707622359309049, 0.0566619690238992, 0.0290912248949643, 0.00957662724216503, 0.0020213587060008, 0.000273561152999412,
        0.000219050647191676, 0.00161857752060636, 0.00766836362659358, 0.0232944318700625, 0.0453713579203497, 0.0566619690238992, 0.0453713579203497, 0.0232944318700625, 0.00766836362659358, 0.00161857752060636, 0.000219050647191676,
        0.000112464352203375, 0.000831005407560634, 0.00393706916085999, 0.0119597601002285, 0.0232944318700625, 0.0290912248949643, 0.0232944318700625, 0.0119597601002285, 0.00393706916085999, 0.000831005407560634, 0.000112464352203375,
        3.70224761237097e-05, 0.000273561152999412, 0.00129605556026987, 0.00393706916085999, 0.00766836362659358, 0.00957662724216503, 0.00766836362659358, 0.00393706916085999, 0.00129605556026987, 0.000273561152999412, 3.70224761237097e-05,
        7.8144113306272e-06, 5.77411237021237e-05, 0.000273561152999412, 0.000831005407560634, 0.00161857752060636, 0.0020213587060008, 0.00161857752060636, 0.000831005407560634, 0.000273561152999412, 5.77411237021237e-05, 7.8144113306272e-06,
        1.05756557075783e-06, 7.8144113306272e-06, 3.70224761237097e-05, 0.000112464352203375, 0.000219050647191676, 0.000273561152999412, 0.000219050647191676, 0.000112464352203375, 3.70224761237097e-05, 7.8144113306272e-06, 1.05756557075783e-06,
        0, 0, 0, 0, 0, 0, 0
    };

    const fp32x8 *weights = (const fp32x8 *) gaussNormalizedWeights;

    const int simdVectorSize = 8;     // 8 floats in the YMM register
    const int windowExtention = 7;    // 7 extra zeroes in the weights array

    const int windowSize = 11;
    const int windowArea = windowSize * windowSize + windowExtention;

    const int chunksPerWindow = windowArea / simdVectorSize;

    // sum has double precision to reduce loss of precision when averaging SSIM over many windows
    double ssimSum = 0;

    // create sliding windows
    fp32x8 window1 [chunksPerWindow];
    fp32x8 window2 [chunksPerWindow];

    // reset last vector in both windows to avoid NaN error in case if there is any garbage
    window1[chunksPerWindow - 1] = window2[chunksPerWindow - 1] = _mm256_setzero_ps();

    // the window moves between [0..(width-11)] and [0..(height-11)]
    for (int window_ypos = 0;  window_ypos <= (imageHeight - windowSize);  window_ypos += 2)
    {
        for (int window_xpos = 0;  window_xpos <= (imageWidth - windowSize);  window_xpos += 2)
        {
            // copy data from image to window arrays

            int internalCounter = 0;
            int currentPixelPosition = window_ypos * imageWidth + window_xpos;

            for (int y = 0;  y < windowSize;  y++)
            {
                for (int x = 0;  x < windowSize;  x++)
                {
                    ((float *) window1)[internalCounter] = (float) image1[currentPixelPosition];
                    ((float *) window2)[internalCounter] = (float) image2[currentPixelPosition];
                    internalCounter++;
                    currentPixelPosition++;
                }
                currentPixelPosition += imageWidth - windowSize;
            }

            // calculate SSIM for a couple of current windows

            // weighted average

            fp32x8 sum1 = _mm256_setzero_ps();
            fp32x8 sum2 = _mm256_setzero_ps();
            for (int i = 0;  i < chunksPerWindow;  i++)
            {
                const register fp32x8 gaussCoefficients = weights[i];
                sum1 += gaussCoefficients * window1[i];
                sum2 += gaussCoefficients * window2[i];
            }

            sum1 += _mm256_permute2f128_ps(sum1, sum1, 1);  // 4 different floats
            sum1 = _mm256_hadd_ps(sum1, sum1);  // 2 different
            const fp32x8 weightedAverage1x8 = _mm256_hadd_ps(sum1, sum1);  // all components contain the same sum

            sum2 += _mm256_permute2f128_ps(sum2, sum2, 1);  // 4 different floats
            sum2 = _mm256_hadd_ps(sum2, sum2);  // 2 different
            const fp32x8 weightedAverage2x8 = _mm256_hadd_ps(sum2, sum2);  // all components contain the same sum

            // weighted variance (squared weighted standard deviation) and weighted covariance

            fp32x8 weightedVariance12x8 = _mm256_setzero_ps();
            fp32x8 weightedCovariancex8 = _mm256_setzero_ps();
            for (int i = 0;  i < chunksPerWindow;  i++)
            {
                const register fp32x8 deviation1 = window1[i] - weightedAverage1x8;
                const register fp32x8 deviation2 = window2[i] - weightedAverage2x8;

                const register fp32x8 gaussCoefficients = weights[i];

                weightedVariance12x8 += gaussCoefficients * (deviation1 * deviation1 + deviation2 * deviation2);
                weightedCovariancex8 += gaussCoefficients * deviation1 * deviation2;
            }

            weightedVariance12x8 += _mm256_permute2f128_ps(weightedVariance12x8, weightedVariance12x8, 1);
            weightedVariance12x8 = _mm256_hadd_ps(weightedVariance12x8, weightedVariance12x8);
            const float weightedVariance12 = _mm_cvtss_f32(_mm256_castps256_ps128( _mm256_hadd_ps(weightedVariance12x8, weightedVariance12x8) ));

            weightedCovariancex8 += _mm256_permute2f128_ps(weightedCovariancex8, weightedCovariancex8, 1);
            weightedCovariancex8 = _mm256_hadd_ps(weightedCovariancex8, weightedCovariancex8);
            const float weightedCovariance = _mm_cvtss_f32(_mm256_castps256_ps128( _mm256_hadd_ps(weightedCovariancex8, weightedCovariancex8) ));

            // SSIM (constants first)

            const float weightedAverage1 = _mm_cvtss_f32(_mm256_castps256_ps128(weightedAverage1x8));
            const float weightedAverage2 = _mm_cvtss_f32(_mm256_castps256_ps128(weightedAverage2x8));

            const float dynamicRange = 255;  // grayscale image, 8 bpp
            const float k1 = 0.01;  // constants from a paper
            const float k2 = 0.03;
            const float c1 = k1 * k1 * dynamicRange * dynamicRange;
            const float c2 = k2 * k2 * dynamicRange * dynamicRange;

            const float ssim = ((2 * weightedAverage1 * weightedAverage2 + c1) * (2 * weightedCovariance + c2)) /
                               ((weightedAverage1 * weightedAverage1 + weightedAverage2 * weightedAverage2 + c1) * (weightedVariance12 + c2));

            // add result to the SSIM accumulator
            ssimSum += ssim;
        }
    }

    // how many 11x11 windows were processed in the image?
    const int numWindows = ((imageHeight - windowSize + 1 + 1) / 2) * ((imageWidth - windowSize + 1 + 1) / 2);

    return ssimSum / numWindows;
}

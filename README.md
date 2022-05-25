# Lua-SigLib
SigLib 

# Frequency Domain Functions
 * SIF_Fft                             - Initialise FFT functionality
 * SDA_Rfft                            - Real to Complex Fast Fourier Transform (FFT)
 * SDA_Cfft                            - Complex to complex FFT
 * SDA_Cifft                           - Complex to complex Inverse FFT
 * SDA_BitReverseReorder               - Bit reverse reorder the data
 * SDA_IndexBitReverseReorder          - Bit reverse reorder the data index
 * SIF_FastBitReverseReorder           - Initialise fast mode bit reverse reordering
 * SDA_RealRealCepstrum                - Real to real cepstrum analysis
 * SDA_RealComplexCepstrum             - Real to complex cepstrum analysis
 * SDA_ComplexComplexCepstrum          - Complex to complex cepstrum analysis
 * SIF_FftTone                         - Initialise FFT tone detection functionality
 * SDA_RfftTone                        - FFT tone detection functionality
 * SDA_Rfftr                           - Real FFT, with real only output
 * SIF_Fft4                            - Initialise radix-4 FFT functionality
 * SDA_Rfft4                           - Real to Complex radix-4 Fast Fourier Transform (FFT)
 * SDA_Cfft4                           - Complex to complex radix-4 FFT
 * SDA_Cifft4                          - Complex to complex Inverse radix-4 FFT
 * SDA_Rfftr4                          - Real radix-4 FFT, with real only output
 * SDA_Cfft2rBy1c                      - Two real FFTs, with one complex FFT
 * SDA_Cfft2rBy1cr                     - Two real FFTs, with one complex FFT, with real only output
 * SDA_Cfft24rBy1c                     - Two real radix-4 FFTs, with one complex radix-4 FFT
 * SDA_Cfft24rBy1cr                    - Two real radix-4 FFTs, with one complex radix-4 FFT, with real only output
 * SDA_BitReverseReorder4              - Bit reverse reorder the radix-4 data
 * SDA_IndexBitReverseReorder4         - Bit reverse reorder the radix-4 data index
 * SIF_FastBitReverseReorder4          - Initialise fast mode radix-4 bit reverse reordering
 * SIF_ZoomFft                         - Initialise zoom FFT functionality
 * SDA_ZoomFft                         - Real to Complex zoom-Fast Fourier Transform (FFT)
 * SIF_ZoomFftSimple                   - Initialise zoom FFT functionality - more efficient but less complex variant
 * SDA_ZoomFftSimple                   - Real to Complex zoom-Fast Fourier Transform (FFT) - more efficient but less complex variant
 * SIF_FdHilbert                       - Initialize frequency domain Hilbert transform function
 * SDA_FdHilbert                       - Frequency domain Hilbert transformer
 * SIF_FdAnalytic                      - Initialize frequency domain analytic signal generator
 * SDA_FdAnalytic                      - Frequency domain analytic signal generator
 * SDA_InstantFreq                     - Instantaneous frequency analyser
 * SDA_Rft                             - Real to Complex Fourier Transform
 * SDA_Rift                            - Real Inverse Fourier Transform
 * SDA_Cft                             - Complex Fourier Transform
 * SDA_Cift                            - Complex Inverse Fourier Transform
 * SDA_FftShift                        - Shift the D.C. bin from 0 to N/2 or V.V. for real data
 * SDA_CfftShift                       - Shift the D.C. bin from 0 to N/2 or V.V. for complex data
 * SDA_FftExtend                       - Extend the array length of a real (or imaginary) FFT array
 * SDA_CfftExtend                      - Extend the array lengths of a real and imaginary FFT arrays
 * SIF_FftArb                          - Initialise the arbitrary FFT operation
 * SUF_FftAllocLengthArb               - Return length of FFT for arbitrary length transform
 * SDA_RfftArb                         - Arbitrary real FFT operation using chirp z-transform
 * SDA_CfftArb                         - Arbitrary complex FFT operation using chirp z-transform
 * SDA_CifftArb                        - Arbitrary complex inverse FFT operation using chirp z-transform
 * SIF_DctII                           - Initialize 1D type II DCT
 * SDA_DctII                           - 1D type II DCT
 * SIF_DctIIOrthogonal                 - Initialize 1D type II orthogonal DCT
 * SDA_DctIIOrthogonal                 - 1D type II orthogonal DCT

# Power Spectrum Functions
 * SIF_FastAutoCrossPowerSpectrum      - Initialise fast auto and cross power spectrum functions
 * SDA_FastAutoPowerSpectrum           - Fast auto power spectrum
 * SDA_FastCrossPowerSpectrum          - Fast cross power spectrum
 * SIF_ArbAutoCrossPowerSpectrum       - Initialise arbitrary length fast auto and cross power spectrum functions
 * SDA_ArbAutoPowerSpectrum            - Arbitrary length auto power spectrum
 * SDA_ArbCrossPowerSpectrum           - Arbitrary length cross power spectrum
 * SIF_WelchPowerSpectrum              - Initilaise the Welch power spectrum functions
 * SDA_WelchRealPowerSpectrum          - Initilaise the Welch power spectrum functions
 * SDA_WelchComplexPowerSpectrum       - Welch complex power spectrum
 * SIF_MagnitudeSquaredCoherence       - Initilaise the magnitude squared coherence functions
 * SDA_MagnitudeSquaredCoherence       - Magnitude squared coherence function

# Frequency Domain Filtering Functions

 * SIF_FirOverlapAdd                   - Initialise overlap and add FFT filtering functionality
 * SDA_FirOverlapAdd                   - Overlap and add FFT filter
 * SIF_FirOverlapSave                  - Initialise overlap and save FFT filtering functionality
 * SDA_FirOverlapSave                  - Overlap and save FFT filter
 * SIF_FftConvolvePre                  - Initialize pre-calculation of frequency domain convolution coefficients
 * SDA_FftConvolvePre                  - Frequency domain convolution with pre-calculated coefficients
 * SDA_FftConvolveArb                  - Frequency domain convolution with arbitrary coefficients
 * SIF_FftCorrelatePre                 - Initialize pre-calculation of frequency domain correlation coefficients
 * SDA_FftCorrelatePre                 - Frequency domain correlation with pre-calculated coefficients
 * SDA_FftCorrelateArb                 - Frequency domain correlation with arbitrary coefficients


# Chirp z-Transform Initialisation Functions

 * SIF_Czt                             - Initialisation function for the chirp z-transform
 * SIF_Awn                             - Generate complex window coeffs
 * SIF_Vl                              - Generate contour definition coeffs
 * SIF_Wm                              - Generate weighting coeffs

# Windowing Functions
 * SIF_Window                          - Initialise windowing functionality
 * SIF_TableTopWindow                  - Initialise the table top (flat centre section) windowing functionality
 * SDA_Window                          - Apply window to a array of data
 * SDA_ComplexWindow                   - Apply a complex window to a complex stream
 * SDA_WindowInverseCoherentGain       - Calculate window function inverse coherent gain
 * SDA_WindowEquivalentNoiseBandwidth  - Calculate window function equivalent noise bandwidth
 * SDA_WindowProcessingGain            - Calculate window function processing gain
 * SDS_I0Bessel                        - Calculate the Izero Bessel function



# Filtering Functions

## FIR Filter Functions

 * SIF_Fir                             - Initialise FIR filter functionality
 * SDS_Fir                             - Perform FIR filter on a data sample
 * SDA_Fir                             - Perform FIR filter on a array of data
 * SDS_FirAddSample                    - Add a sample to the filter delay line
 * SDA_FirAddSamples                   - Add multiple samples to the filter delay line
 * SIF_Comb                            - Initialise comb filter functionality
 * SDS_Comb                            - Perform a comb filter on a data sample
 * SDA_Comb                            - Perform a comb filter on a data array
 * SIF_FirComplex                      - Initialise a complex filter
 * SDS_FirComplex                      - Perform a complex filter on a data sample
 * SDA_FirComplex                      - Perform a complex filter on a data stream
 * SIF_FirWithStore                    - Initialise FIR filter functionality - With store maybe more efficient than standard FIR filter functions
 * SDS_FirWithStore                    - Perform FIR filter on a data sample - With store maybe more efficient than standard FIR filter functions
 * SDA_FirWithStore                    - Perform FIR filter on a array of data - With store maybe more efficient than standard FIR filter functions
 * SIF_FirComplexWithStore             - Initialise complex FIR filter functionality - With store maybe more efficient than standard FIR filter functions
 * SDS_FirComplexWithStore             - Perform complex FIR filter on a data sample - With store maybe more efficient than standard FIR filter functions
 * SDA_FirComplexWithStore             - Perform complex FIR filter on a array of data - With store maybe more efficient than standard FIR filter functions
 * SDS_FirWithStoreAddSample           - Add a sample to the filter with store delay line
 * SDA_FirWithStoreAddSamples          - Add multiple samples to the filter with store delay line
 * SIF_FirExtendedArray                - Initialise FIR filter functionality - With extended state array maybe more efficient than standard FIR filter functions
 * SDS_FirExtendedArray                - Perform FIR filter on a data sample - With extended state array maybe more efficient than standard FIR filter functions
 * SDA_FirExtendedArray                - Perform FIR filter on a array of data - With extended state array maybe more efficient than standard FIR filter functions
 * SIF_FirComplexExtendedArray         - Initialise a complex filter - With extended state array maybe more efficient than standard FIR filter functions
 * SDS_FirComplexExtendedArray         - Perform a complex filter on a data sample - With extended state array maybe more efficient than standard FIR filter functions
 * SDA_FirComplexExtendedArray         - Perform a complex filter on a data stream - With extended state array maybe more efficient than standard FIR filter functions
 * SDS_FirExtendedArrayAddSample       - Add a sample to the filter with extended state array
 * SDA_FirExtendedArrayAddSamples      - Add multiple samples to the filter with extended state array
 * SIF_LowPassFirFilter                - Create a low pass filter with the given centre frequency
 * SIF_HighPassFirFilter               - Create a high pass filter with the given centre frequency
 * SIF_BandPassFirFilter               - Create a band pass filter with the given centre frequency
 * SIF_LowPassFirFilterWindow          - Create a low pass filter with the given centre frequency, using windowing
 * SIF_HighPassFirFilterWindow         - Create a high pass filter with the given centre frequency, using windowing
 * SIF_BandPassFirFilterWindow         - Create a band pass filter with the given centre frequency, using windowing
 * SUF_KaiserApproximation             - Approximate the desired number of FIR filter coefficients using the Kaiser approximation
 * SIF_MatchedFirFilter                - Initialize coefficients for a matched FIR filter from a given data set
 * SDA_FirFilterInverseCoherentGain    - Calculate the inverse coherent gain for an FIR filter
 * SIF_TappedDelayLine                 - Initialize the tapped delay line / multi-path functions
 * SDS_TappedDelayLine                 - Tapped delay line / multi-path function on a per sample basis
 * SDA_TappedDelayLine                 - Tapped delay line / multi-path function on a per array basis
 * SIF_TappedDelayLineComplex          - Initialize the complex tapped delay line / multi-path functions
 * SDS_TappedDelayLineComplex          - Complex tapped delay line / multi-path function on a per sample basis
 * SDA_TappedDelayLineComplex          - Complex tapped delay line / multi-path function on a per array basis
 * SIF_TappedDelayLineIQ               - Initialize the IQ tapped delay line / multi-path functions
 * SDS_TappedDelayLineIQ               - IQ tapped delay line / multi-path function on a per sample basis
 * SDA_TappedDelayLineIQ               - IQ tapped delay line / multi-path function on a per array basis
 * SIF_FirPolyPhaseGenerate            - Initialize a polyphase FIR filter

## IIR Filter Functions

 * SIF_Iir                             - Initialise IIR filter functionality - uses the biquad filter structure - subtraction of feedback coefficients
 * SDS_Iir                             - Perform IIR filter on a data sample - uses the biquad filter structure - subtraction of feedback coefficients
 * SDA_Iir                             - Perform IIR filter on a array of data - uses the biquad filter structure - addition of negated feedback coefficients
 * SDS_IirMac                          - Perform IIR filter on a data sample - uses the biquad filter structure - addition of negated feedback coefficients
 * SDA_IirMac                          - Perform IIR filter on a array of data - uses the biquad filter structure
 * SIF_IirOrderN                       - Initialise the order N IIR filter functionality - uses a single order N structure
 * SDS_IirOrderN                       - Perform an order N IIR filter on a data sample - uses a single order N structure
 * SDA_IirOrderN                       - Perform an order N IIR filter on a array of data - uses a single order N structure
 * SIF_NcIir                           - Initialise bi-directional (non-causal) IIR filter functionality
 * SDA_NcIir                           - Perform bi-directional (non-causal) IIR filter on a array of data
 * SDA_BilinearTransform               - Transform poles and zeros from s-plane to z-plane using the bilinear transform
 * SDA_PreWarp                         - Pre-warps frequencies for bilinear transform
 * SDA_MatchedZTransform               - Transform poles and zeros from s-plane to z-plane using the matched z-transform
 * SDA_IirZplaneToCoeffs               - Translate rectangular z-plane poles and zeros to IIR coefficients
 * SDA_IirZplanePolarToCoeffs          - Translate polar z-plane poles and zeros to IIR coefficients
 * SDA_IirZplaneLpfToLpf               - Frequency translate a low-pass IIR filter
 * SDA_IirZplaneLpfToHpf               - Translate a low-pass IIR filter to high pass
 * SDA_IirZplaneLpfToBpf               - Translate a low-pass IIR filter to band pass
 * SDA_IirZplaneLpfToBsf               - Translate a low-pass IIR filter to band stop
 * SDA_ModifyIirFilterGain             - Modify the pass band gain of the IIR filter
 * SIF_LowPassIirFilter                - Initialize the coefficients of a low-pass IIR filer biquad
 * SIF_HighPassIirFilter               - Initialize the coefficients of a high-pass IIR filer biquad
 * SIF_AllPassIirFilter                - Initialize the coefficients of an all-pass IIR filer biquad
 * SIF_BandPassIirFilter               - Initialize the coefficients of a band-pass IIR filer biquad
 * SIF_NotchIirFilter                  - Initialize the coefficients of a notch IIR filer biquad
 * SIF_PeakingIirFilter                - Initialize the coefficients of a peaking IIR filer biquad
 * SIF_LowShelfIirFilter               - Initialize the coefficients of a low-shelf IIR filer biquad
 * SIF_HighShelfIirFilter              - Initialize the coefficients of a high-shelf IIR filer biquad
 * SDS_RemoveDC                        - Remove the DC component of a signal using an IIR filter
 * SDA_RemoveDC                        - Remove the DC component of a signal using an IIR filter
 * SIF_OnePole                         - Initialize the one pole low-pass filter functions
 * SDS_OnePole                         - Apply a one pole low-pass filter on a data sample
 * SDA_OnePole                         - Apply a one pole low-pass filter on an array of data
 * SDS_OnePoleNormalized               - Apply a one pole low-pass filter on a data sample. The step response gain is normalized to 1.0
 * SDA_OnePoleNormalized               - Apply a one pole low-pass filter on an array of data. The step response gain is normalized to 1.0
 * SDA_OnePolePerSample                - Apply a one pole low-pass filter, between arrays of data
 * SIF_OnePoleHighPass                 - Initialize one pole high-pass filter functions
 * SDS_OnePoleHighPass                 - Apply a one pole high-pass filter on a data sample
 * SDA_OnePoleHighPass                 - Apply a one pole high-pass filter on an array of data
 * SDS_OnePoleHighPassNormalized       - Apply a one pole high-pass filter on a data sample. The step response gain is normalized to 1.0
 * SDA_OnePoleHighPassNormalized       - Apply a one pole high-pass filter on an array of data. The step response gain is normalized to 1.0
 * SDA_OnePoleHighPassPerSample        - Apply a one pole high-pass filter, between arrays of data
 * SDS_OnePoleHighPassCutOffFrequencyToFilterCoeff - Convert high-pass cut-off frequency to one-pole coefficient
 * SDS_TimeConstantToOnePoleFilterCoeff    - Convert a time constant to a one-pole IIR filter coefficient
 * SDS_CutOffFrequencyToOnePoleFilterCoeff - Convert a cut-off frequency to a one-pole IIR filter coefficient
 * SIF_AllPole                         - Initialise the all pole IIR filter functionality
 * SDS_AllPole                         - Apply an all pole filter on a data sample
 * SDA_AllPole                         - Apply an all pole filter on an array of data
 * SDA_ZDomainCoefficientReorg         - Reorganises z-domain coefficients from Digital Filter Plus
 * SIF_IirNotchFilter                  - Generate the coefficients for an IIR notch filter
 * SIF_IirNormalizedCoefficients       - Generate normalized Butterworth and Bessel low-pass filter coefficients
 * SIF_IirNormalizedSPlaneCoefficients - Generate normalized s-Plane Butterworth and Bessel low-pass filter coefficients
 * SDA_TranslateSPlaneCutOffFrequency  - Translate the cut-off frequency of a low-pass filter
 * SDA_IirLpLpShift                    - Shift the cut-off frequency of an IIR digital filter
 * SDA_IirLpHpShift                    - Transform a LPF IIR filter into a HPF and shift the cut-off frequency
 * SIF_Iir2PoleLpf                     - Two pole IIR low-pass filter design
 * SDS_Iir2Pole                        - Apply two pole IIR filter on a data sample
 * SDA_Iir2Pole                        - Apply two pole IIR filter on an array of data
 * SDA_IirNegateAlphaCoeffs            - Negate IIR filter feedback coefficients to use the MAC variants - these may be more efficient on certain processors

## Generic Filtering Functions
 * SDA_Integrate                       - Integrate the waveform
 * SDA_Differentiate                   - Differentiate the waveform
 * SIF_LeakyIntegrator                 - Initialise the leaky integrator
 * SDS_LeakyIntegrator1                - Leaky integrator - Limit and add input
 * SDS_LeakyIntegrator2                - Leaky integrator - Add input and limit
 * SIF_HilbertTransformer              - Initialise Hilbert transformer filter coefficients
 * SIF_GoertzelFilter                  - Initialise Goertzel filter coefficients
 * SDA_GoertzelFilter                  - Apply Goertzel filter
 * SDS_GoertzelFilter                  - Apply Goertzel filter (per sample)
 * SIF_GoertzelDetect                  - Initialise Goertzel filter in signal detect mode
 * SDA_GoertzelDetect                  - Apply Goertzel filter in signal detect mode
 * SIF_GoertzelFilterComplex           - Initialise Goertzel filter coefficients
 * SDA_GoertzelDetectComplex           - Apply Goertzel filter in signal detect mode
 * SIF_GaussianFilter                  - Initialise Gaussian filter coefficients
 * SIF_GaussianFilter2                 - Initialise Gaussian filter coefficients - second type
 * SIF_RaisedCosineFilter              - Initialise a raised cosine filter coefficients
 * SIF_RootRaisedCosineFilter          - Initialise a root raised cosine filter coefficients
 * SDS_ZTransform                      - Calculate the z-Transform of the source array
 * SDS_ZTransformDB                    - Calculate the z-Transform of the source array in dB
 * SUF_EstimateBPFilterLength          - Estimate the length of a band-pass FIR filer
 * SUF_EstimateBPFilterError           - Estimate the error of a given band-pass FIR filter
 * SUF_FrequenciesToOctaves            - Convert frequencies to octaves for given sample rate
 * SUF_FrequenciesToCentreFreqHz       - Convert frequencies to centre frequencies for given sample rate
 * SUF_FrequenciesToQFactor            - Convert frequencies to Q Factor
 * SUF_BandwidthToQFactor              - Convert bandwidth to Q factor
 * SUF_QFactorToBandwidth              - Convert Q factor to bandwidth

# Acoustic Processing Functions

 * SDA_LinearMicrophoneArrayBeamPattern        - Beam pattern for a linear microphone array
 * SDA_LinearMicrophoneArrayBeamPatternLinear  - Beam pattern for a linear microphone array
 * SDA_MicrophoneArrayBeamPattern      - Beam pattern for a 2D microphone array
 * SDA_MicrophoneArrayBeamPatternLinear- Beam pattern for a 2D microphone array
 * SDA_MicrophoneArrayCalculateDelays  - Calculate delays for direction for a microphone array
 * SDA_LinearMicrophoneArrayBeamPattern- Calculate the beam pattern for a linear microphone array
 * SDA_LinearMicrophoneArrayBeamPattern- Calculate the beam pattern for a linear microphone array
 * SDA_MicrophoneArrayBeamPattern      - Calculate the beam pattern for an arbitrary microphone array
 * SDS_TemperatureToSpeedOfSoundInAir  - Adjust speed of sound in air for temperature variation


# Adaptive Filtering Functions


 * SIF_Lms                             - Initialise adaptive filter functionality
 * SDS_Lms                             - Perform filter on a data sample
 * SDA_LmsUpdate                       - Update LMS filter taps
 * SDA_LeakyLmsUpdate                  - Update leaky LMS filter taps
 * SDA_NormalizedLmsUpdate             - Update normalized LMS filter taps
 * SDA_SignErrorLmsUpdate              - Update signed error LMS filter taps
 * SDA_SignDataLmsUpdate               - Update signed data LMS filter taps
 * SDA_SignSignLmsUpdate               - Update signed sign LMS filter taps

# Convolution, Deconvolution And Correlation Functions


 * SDA_ConvolveLinear                  - Linearly convolve input with impulse response
 * SDA_ConvolvePartial                 - Partially convolve input with impulse response
 * SDA_ConvolveCircular                - Cyclic convolve input with impulse response
 * SDA_ConvolveLinearComplex           - Linearly convolve complex input with complex impulse response
 * SDA_ConvolvePartialComplex          - Partially convolve complex input with complex impulse response
 * SDA_ConvolveCircularComplex         - Cyclic convolve complex input with complex impulse response
 * SDA_FftDeconvolution                - Deconvolve a signal and an impulse response
 * SIF_FftDeconvolutionPre             - Initialise the FFT deconvolution function, with pre-calculated FFT of impulse response
 * SDA_FftDeconvolutionPre             - Deconvolve a signal and an impulse response, with pre-calculated FFT of impulse response
 * SDA_CorrelateLinear                 - Linear cross correlate two data arrays
 * SDA_CorrelatePartial                - Partial linear cross correlate two data arrays
 * SDA_CorrelateCircular               - Cyclic cross correlate two data arrays
 * SDA_Covariance                      - Calculate the covariance of two data arrays
 * SDA_CovariancePartial               - Calculate part of the covariance of two data arrays
 * SDA_CorrelateLinearReturnPeak       - Linear cross correlate two data arrays and return the peak


# Delay Functions

 * SIF_FixedDelay                      - Initialise the fixed delay functions
 * SDS_FixedDelay                      - Delay the data by a fixed number of samples
 * SDA_FixedDelay                      - Delay the data in an array by a fixed number of samples
 * SIF_FixedDelayComplex               - Initialise the fixed complex delay functions
 * SDS_FixedDelayComplex               - Delay the complex data by a fixed number of samples
 * SDA_FixedDelayComplex               - Delay the complex data in an array by a fixed number of samples
 * SDA_ShortFixedDelay                 - Delay the data in an array by a fixed small number of samples
 * SIF_VariableDelay                   - Initialise the variable delay functions
 * SDS_VariableDelay                   - Delay the data by a variable number of samples
 * SDA_VariableDelay                   - Delay the data in an array by a variable number of samples
 * SIF_VariableDelayComplex            - Initialise the variable complex delay functions
 * SDS_VariableDelayComplex            - Delay the complex data by a variable number of samples
 * SDA_VariableDelayComplex            - Delay the complex data in an array by a variable number of samples
 * SUF_IncreaseVariableDelay           - Increase the variable delay by a number of samples
 * SUF_DecreaseVariableDelay           - Decrease the variable delay by a number of samples
 * SDA_Align                           - Align two arrays of data


# Image Processing Functions

 * SIM_Fft2d                           - Perform 2d FFT on image
 * SIF_Fft2d                           - Initialise 2d FFT
 * SIM_Conv3x3                         - Convolve an image with a 3x3 kernel
 * SIM_Sobel3x3                        - Apply a 3x3 Sobel edge detection filter
 * SIM_SobelVertical3x3                - Apply a vertical 3x3 Sobel filter
 * SIM_SobelHorizontal3x3              - Apply a horizontal 3x3 Sobel filter
 * SIM_Median3x3                       - Apply a 3x3 median filter
 * SIF_ConvCoefficients3x3             - Initialize a 3x3 convolution filter
 * SIM_Max                             - Return the maximum value in the image
 * SIM_Min                             - Return the minimum value in the image
 * SDA_Histogram                       - Produce histogram of data
 * SDA_HistogramEqualize               - Equalize the histogram of the data


# Machine Learning Functions


 * SDA_TwoLayer2CategoryNetworkFit     - Two layer, two category network fit (train)
 * SDA_TwoLayer2CategoryNetworkPredict - Two layer, two category network predict (infer)
 * SDA_TwoLayerNCategoryNetworkFit     - Two layer, N category network fit (train)
 * SDA_TwoLayerNCategoryNetworkPredict - Two layer, N category network predict (infer)
 * SDS_ActivationReLU                  - ReLU activation function (sample oriented)
 * SDA_ActivationReLU                  - ReLU activation function (array oriented)
 * SDS_ActivationReLUDerivative        - ReLU derivative activation function (sample oriented)
 * SDA_ActivationReLUDerivative        - ReLU derivative activation function (array oriented)
 * SDS_ActivationLeakyReLU             - Leaky ReLU activation function (sample oriented)
 * SDA_ActivationLeakyReLU             - Leaky ReLU activation function (array oriented)
 * SDS_ActivationLeakyReLUDerivative   - Leaky ReLU derivative activation function (sample oriented)
 * SDA_ActivationLeakyReLUDerivative   - Leaky ReLU derivative activation function (array oriented)
 * SDS_ActivationLogistic              - Logistic activation function (sample oriented)
 * SDA_ActivationLogistic              - Logistic activation function (array oriented)
 * SDS_ActivationLogisticDerivative    - Logistic derivative activation function (sample oriented)
 * SDA_ActivationLogisticDerivative    - Logistic derivative activation function (array oriented)
 * SDS_ActivationTanH                  - Hyperbolic tangent activation function (sample oriented)
 * SDA_ActivationTanH                  - Hyperbolic tangent activation function (array oriented)
 * SDS_ActivationTanHDerivative        - Hyperbolic tangent derivative activation function (sample oriented)
 * SDA_ActivationTanHDerivative        - Hyperbolic tangent derivative activation function (array oriented)

# Image Coding Functions

 * SIF_Dct8x8                          - Initialise the DCT functionality
 * SIM_Dct8x8                          - Apply an 8x8 DCT to an image block
 * SIM_Idct8x8                         - Apply an inverse 8x8 DCT to an image block
 * SIM_ZigZagScan                      - Zig-zag scan an image block
 * SIM_ZigZagDescan                    - Zig-zag scan an image array


# Signal Generation Functions


 * SDA_SignalGenerate                  - Generate an array full of a signal (sin, cos, square, ramp etc.)
 * SDS_SignalGenerate                  - Generate a sample of a signal (sin, cos, square, ramp etc.)
 * SDA_SigGenRamp                      - Generate a ramp signal
 * SIF_Resonator                       - Initialise the resonator
 * SDA_Resonator                       - A digital resonator
 * SIF_Resonator                       - Initialise the resonator1
 * SDA_Resonator1                      - Second digital resonator
 * SDA_Resonator1Add                   - Add a sample into digital resonator1
 * SDA_SignalGeneratePolarWhiteNoise   - Add a polar white noise signal to the input array
 * SDS_SignalGeneratePolarWhiteNoise   - Add a polar white noise signal to the input sample
 * SDA_SignalGeneratePolarGaussianNoise - Add a polar Gaussian white noise signal to the input array
 * SDS_SignalGeneratePolarGaussianNoise - Add a polar Gaussian white noise signal to the input sample
 * SDA_SignalAddPolarJitterAndGaussianNoise - Add jitter and a Gaussian polar white noise signal to the input array
 * SDS_SignalAddPolarJitterAndGaussianNoise - Add jitter and a Gaussian polar white noise signal to the input sample
 * SDA_Ramp                              - Generate a ramp in the array
 * SIF_RandomNumber                      - Initialize random number seed
 * SDS_RandomNumber                      - Generate a random number
 * SDA_RandomNumber                      - Fill the array with random numbers


# Modulation And Communication Functions


 * SDA_BitErrorRate                    - Calculate bit error rate for a signal
 * SDA_Interleave                      - Interleave the samples in a data stream
 * SDA_Deinterleave                    - De-interleave the samples in a data stream
 * SCV_EuclideanDistance               - Return the Euclidean distance between two vectors
 * SCV_EuclideanDistanceSquared        - Return the square of the Euclidean distance between two vectors
 * SDS_EuclideanDistance               - Return the Euclidean distance between two samples
 * SDS_EuclideanDistanceSquared        - Return the square of the Euclidean distance between two samples
 * SDA_EuclideanDistance               - Return the Euclidean distance between two arrays of samples
 * SDA_EuclideanDistanceSquared        - Return the square of the Euclidean distance between two arrays of samples
 * SCA_EuclideanDistance               - Return the Euclidean distance between two arrays of vectors
 * SCA_EuclideanDistanceSquared        - Return the square of the Euclidean distance between two arrays of vectors
 * SDS_ManchesterEncode                - Apply Manchester encoding to the source data bit
 * SDS_ManchesterDecode                - Apply Manchester decoding to the source data bit
 * SDS_ManchesterEncodeByte            - Apply Manchester encoding to the source data byte
 * SDS_ManchesterDecodeByte            - Apply Manchester decoding to the source data byte

 * SIF_DetectNumericalWordSequence     - Initialize function to detect a sequence of numerical words
 * SDS_DetectNumericalWordSequence     - Detect a sequence of numerical words
 * SIF_DetectNumericalBitSequence      - Initialize function to detect a sequence of numerical bits
 * SDS_DetectNumericalBitSequence      - Detect a sequence of numerical bits
 * SIF_DetectCharacterSequence         - Initialize function to detect a sequence of characters
 * SDS_DetectCharacterSequence         - Detect a sequence of characters
 * SDS_ErrorVector                     - Calculate error vector on the input samples
 * SDS_ErrorVectorMagnitudePercent     - Calculate error vector magnitude in percentage on the input samples
 * SDS_ErrorVectorMagnitudeDecibels    - Calculate error vector magnitude in dB on the input samples
 * SDS_ReverseDiBits                   - Reverse the di-bit pair to swap from "computer" order to "ITU-T" order
 * SDS_QpskBitErrorCount               - Count the bit errors in the QPSK data stream
 * SDS_BitErrorRate                    - Calculate the bit error rate in the QPSK data stream

 * SIF_PhaseLockedLoop                 - Initialise phase locked loop
 * SDS_PhaseLockedLoop                 - Phase locked loop on per sample basis
 * SDA_PhaseLockedLoop                 - Phase locked loop on array basis
 * SIF_CostasLoop                      - Initialise Costas loop
 * SDS_CostasLoop                      - Costas loop on per sample basis
 * SDA_CostasLoop                      - Costas loop on array basis
 * SRF_CostasLoop                      - Reset the Costas loop functions
 * SIF_180DegreePhaseDetect            - Initialise the 180 degree phase shift detector
 * SDA_180DegreePhaseDetect            - Detect 180 degree phase shifts in the input signal

 * SIF_TriggerReverberator             - Initialise the trigger reverberator - ensures trigger continues if signal diminishes
 * SDA_TriggerReverberator             - Apply the trigger reverberator on per sample basis
 * SDS_TriggerReverberator             - Apply the trigger reverberator on array basis
 * SDA_TriggerSelector                 - Trigger selector function

 * SIF_EarlyLateGate                   - Initialise the early-late-gate timing error detector
 * SDA_EarlyLateGate                   - Apply the early-late-gate timing error detector on per sample basis
 * SDS_EarlyLateGate                   - Apply the early-late-gate timing error detector on array basis
 * SIF_EarlyLateGateSquarePulse        - Initialise the early-late-gate timing error detector on a square pulse input
 * SDA_EarlyLateGateSquarePulse        - Apply the early-late-gate timing error detector on per sample basis
 * SDS_EarlyLateGateSquarePulse        - Apply the early-late-gate timing error detector on array basis

 * SDS_ConvEncoderK3                   - K = 3, rate 1/2 convolutional encoder (trellis coded modulator)
 * SIF_ViterbiDecoderK3                - Initialise K = 3, rate 1/2 Viterbi decoder
 * SDS_ViterbiDecoderK3                - K = 3, rate 1/2 Viterbi decoder
 * SDS_ConvEncoderV32                  - V.32 (32QAM) convolutional encoder (trellis coded modulator)
 * SIF_ViterbiDecoderV32               - Initialise V.32 (32QAM) Viterbi decoder
 * SDS_ViterbiDecoderV32               - V.32 (32QAM) Viterbi decoder

 * SIF_AmplitudeModulate               - Initialize amplitude modulation function
 * SDA_AmplitudeModulate               - Amplitude modulate a signal on a per array basis
 * SDS_AmplitudeModulate               - Amplitude modulate a signal on a per sample basis
 * SIF_AmplitudeModulate2              - Initialize amplitude modulation function - version 2
 * SDA_AmplitudeModulate2              - Amplitude modulate a signal on a per array basis - version 2
 * SDS_AmplitudeModulate2              - Amplitude modulate a signal on a per sample basis - version 2
 * SIF_ComplexShift                    - Initialise complex frequency shifting
 * SDA_ComplexShift                    - Complex frequency shift a signal. Can also operate as a lock-in amplifier
 * SIF_FrequencyModulate               - Initialize frequency modulation functions
 * SDS_FrequencyModulate               - Frequency modulate a signal - per sample basis
 * SDA_FrequencyModulate               - Frequency modulate a signal - array basis
 * - Can also be used for Voltage Controlled Oscillator (VCO or NCO)
 * SDA_FrequencyDemodulate             - Demodulate an FM signal
 * SIF_FrequencyModulateComplex        - Initialize complex frequency modulation functions
 * SDS_FrequencyModulateComplex        - Frequency modulate a complex signal - per sample basis
 * SDA_FrequencyModulateComplex        - Frequency modulate a complex signal - array basis
 * SDA_DeltaModulate                   - Delta modulate a signal
 * SDA_DeltaDemodulate                 - Demodulate an delta modulated signal
 * SDA_DeltaModulate2                  - Another function to delta modulate a signal

 * SIF_CostasQamDemodulate             - Initilize the Costas loop based QAM / QPSK demodulator functions. Includes early-late-gate timing error detector
 * SDS_CostasQamDemodulate             - Implement the Costas loop based QAM / QPSK demodulator functions on a per sample basis
 * SDS_CostasQamDemodulateDebug        - Implement the Costas loop based QAM / QPSK demodulator functions on a per sample basis with debug information
 * SDA_CostasQamDemodulate             - Implement the Costas loop based QAM / QPSK demodulator functions on an array basis
 * SDA_CostasQamDemodulateDebug        - Implement the Costas loop based QAM / QPSK demodulator functions on an array basis with debug information
 * SIF_QpskModulate                    - Initialise the QPSK modulation function
 * SDA_QpskModulate                    - QPSK modulate a signal
 * SIF_QpskDemodulate                  - Initialise the QPSK demodulation function
 * SDA_QpskDemodulate                  - Demodulate a QPSK signal
 * SDA_QpskDemodulateDebug             - Demodulate a QPSK signal, with debug output
 * SDA_QpskDifferentialEncode          - QPSK differential encoding
 * SDA_QpskDifferentialDecode          - QPSK differential decoding
 * SIF_FskModulate                     - Initialise the FSK modulation function
 * SDA_FskModulateByte                 - FSK modulate a signal - 8 bit byte input
 * SDA_FskDemodulateByte               - Demodulate a FSK (and CP-FSK) signal - 8 bit byte output
 * SDA_CpfskModulateByte               - Continuous phase FSK modulate a signal - 8 bit byte input
 * SDA_FskModulate                     - FSK modulate a signal - 1 bit input
 * SDA_FskDemodulate                   - Demodulate an FSK (and CP-FSK) signal - 1 bit output
 * SDA_CpfskModulate                   - Continuous phase FSK modulate a signal - 1 bit input
 * SIF_Qam16Modulate                   - Initialise the QAM-16 modulation function
 * SDA_Qam16Modulate                   - QAM-16 modulate a signal
 * SIF_Qam16Demodulate                 - Initialise the QAM-16 demodulation function
 * SDA_Qam16Demodulate                 - QAM-16 demodulate a signal
 * SDA_Qam16DemodulateDebug            - QAM-16 demodulate a signal (with debug information)
 * SDA_Qam16DifferentialEncode         - QAM-16 differential encoding
 * SDA_Qam16DifferentialDecode         - QAM-16 differential decoding
 * SIF_BpskModulate                    - Initialise the BPSK modulation function
 * SDA_BpskModulate                    - BPSK modulate a signal
 * SDA_BpskModulateByte                - BPSK modulate a signal - 8 bit byte input
 * SIF_BpskDemodulate                  - Initialise the BPSK demodulation function
 * SDA_BpskDemodulate                  - Demodulate a BPSK signal
 * SDA_BpskDemodulateDebug             - Demodulate a BPSK signal, with debug output
 * SIF_DpskModulate                    - Initialise the DPSK modulation function
 * SDA_DpskModulate                    - DPSK modulate a signal
 * SDA_DpskModulateByte                - DPSK modulate a signal - 8 bit byte input
 * SIF_DpskDemodulate                  - Initialise the DPSK demodulation function
 * SDA_DpskDemodulate                  - Demodulate a DPSK signal
 * SDA_DpskDemodulateDebug             - Demodulate a DPSK signal, with debug output
 * SIF_PiByFourDQpskModulate           - Initilize the PI/4 D-QPSK modulation function
 * SDA_PiByFourDQpskModulate           - PI/4 D-QPSK modulate a signal
 * SDS_ChannelizationCode              - 3GPP2 compliant channelization code generation
 * SDA_ComplexQPSKSpread               - 3GPP2 compliant complex weighting, spreading and scrambling
 * SDA_ComplexQPSKDeSpread             - 3GPP2 compliant complex de-weighting, de-spreading and de-scrambling

 * SUF_AsyncCharacterLength            - Return the asynchronous character length for a given word length and number of start, stop and partiy bits
 * SDA_SyncToAsyncConverter            - Convert a data sequence from synchronous to asynchronous
 * SDA_AsyncToSyncConverter            - Convert a data sequence from asynchronous to synchronous
 * SIF_AsyncAddRemoveStopBits          - Initilize the function to remove stop bits from an asynchronous sequence
 * SDA_AsyncRemoveStopBits             - Remove stop bits from an asynchronous sequence
 * SDA_AsyncAddStopBits                - Add stop bits to an asynchronous sequence
 * SDA_DecreaseWordLength              - Decrease the wordlength of a synchronous sequence
 * SDA_IncreaseWordLength              - Increase the wordlength of a synchronous sequence

 * SDS_Scrambler1417                   - 1 + x-14 + x-17 polynomial Pseudo Random Binary Sequence (PRBS) Cyclic Redundancy Check (CRC) scrambler
 * SDS_Descrambler1417                 - 1 + x-14 + x-17 polynomial Pseudo Random Binary Sequence (PRBS) Cyclic Redundancy Check (CRC) descrambler
 * SDS_Scrambler1417WithInversion      - 1 + x-14 + x-17 polynomial Pseudo Random Binary Sequence (PRBS) Cyclic Redundancy Check (CRC) scrambler - with inversion
 * SDS_Descrambler1417WithInversion    - 1 + x-14 + x-17 polynomial Pseudo Random Binary Sequence (PRBS) Cyclic Redundancy Check (CRC) descrambler - with inversion
 * SDS_Scrambler1823                   - 1 + x-18 + x-23 polynomial Pseudo Random Binary Sequence (PRBS) Cyclic Redundancy Check (CRC) scrambler
 * SDS_Descrambler1823                 - 1 + x-18 + x-23 polynomial Pseudo Random Binary Sequence (PRBS) Cyclic Redundancy Check (CRC) descrambler
 * SDS_Scrambler523                    - 1 + x-5 + x-23 polynomial Pseudo Random Binary Sequence (PRBS) Cyclic Redundancy Check (CRC) scrambler
 * SDS_Descrambler523                  - 1 + x-5 + x-23 polynomial Pseudo Random Binary Sequence (PRBS) Cyclic Redundancy Check (CRC) descrambler
 * SDS_ScramblerDescramblerPN9         - PN-9 scrambler / descrambler
 * SDS_ScramblerDescramblerPN15        - PN-15 scrambler / descrambler
 * SDS_SequenceGeneratorPN9            - PN-9 scrambler / descrambler sequence generator
 * SDS_SequenceGeneratorPN15           - PN-9 scrambler / descrambler sequence generator
 * SDS_ScramblerDescramblergCRC24      - 3GPP CRC24 descrambler
 * SDS_SequenceGeneratorgCRC24         - 3GPP CRC24 sequence generator
 * SDS_ScramblerDescramblergCRC16      - 3GPP CRC16 descrambler
 * SDS_SequenceGeneratorgCRC16         - 3GPP CRC16 sequence generator
 * SDS_ScramblerDescramblergCRC12      - 3GPP CRC12 descrambler
 * SDS_SequenceGeneratorgCRC12         - 3GPP CRC12 sequence generator
 * SDS_ScramblerDescramblergCRC8       - 3GPP CRC8 descrambler
 * SDS_SequenceGeneratorgCRC8          - 3GPP CRC8 sequence generator
 * SDS_LongCodeGenerator3GPPDL         - 3GPP downlink long code sequence generator
 * SDS_LongCodeGenerator3GPPUL         - 3GPP uplink long code sequence generator
 * SDA_Multiplex                       - Multiplex the elements in the arrays
 * SDA_Demultiplex                     - Demultiplex the elements in the arrays
 * SDA_Mux_N                           - Multiplex N channels onto one stream
 * SDA_Demux_N                         - Demultiplex N channels from one steam


# Interpolation And Decimation Functions


 * SIF_Decimate                        - Initialise the decimation function
 * SDA_Decimate                        - Decimate the sample rate of the source data
 * SIF_Interpolate                     - Initialise the interpolation function
 * SDA_Interpolate                     - Interpolate the sample rate of the source data
 * SIF_FilterAndDecimate               - Initialise the filter and decimation function
 * SDA_FilterAndDecimate               - Filter and decimate the sample rate of the source data
 * SIF_InterpolateAndFilter            - Initialise the interpolate and filter function
 * SDA_InterpolateAndFilter            - Interpolate and filter sample rate of the source data
 * SDA_ResampleLinear                  - Linear re-sampler
 * SDA_ResampleLinearNSamples          - Linear re-sampler and output N samples
 * SDA_InterpolateLinear1D             - Linear interpolation version 1D
 * SDA_InterpolateLinear2D             - Linear interpolation version 2D
 * SIF_ResampleSinc                    - Initialize sin(x)/x (sinc) re-sampler
 * SIF_ResampleWindowedSinc            - Initialize window'd sin(x)/x (sinc) re-sampler
 * SDA_ResampleSinc                    - sin(x)/x (sinc) re-sampler
 * SDA_ResampleSincNSamples            - sin(x)/x (sinc) re-sampler and output N samples
 * SIF_InterpolateSinc1D               - Initialize sin(x)/x (sinc) interpolation
 * SIF_InterpolateWindowedSinc1D       - Initialize window'd sin(x)/x (sinc) interpolation
 * SDA_InterpolateSinc1D               - sin(x)/x (sinc) interpolation
 * SIF_ResampleLinearContiguous        - Initialize contiguous (over entire array) linear re-sampler
 * SDA_ResampleLinearContiguous        - Contiguous (over entire array) linear re-sampler
 * SIF_ResampleSincContiguous          - Initialize contiguous (over entire array) sin(x)/x (sinc) re-sampler
 * SIF_ResampleWindowedSincContiguous  - Initialize contiguous (over entire array) sin(x)/x (sinc) re-sampler
 * SDA_ResampleSincContiguous          - Contiguous (over entire array) sin(x)/x (sinc) re-sampler

 * SDS_InterpolateQuadratic1D          - Quadratic interpolate the sample interpolation
 * SDS_InterpolateQuadraticBSpline1D   - B-Spline interpolate the sample interpolation
 * SDS_InterpolateQuadraticLagrange1D  - Lagrange interpolate the sample interpolation

# DTMF Generation And Detection Functions


 * SIF_DtmfGenerate                    - Initialise the DTMF generation function
 * SDA_DtmfGenerate                    - Generate the DTMF signals
 * SIF_DtmfDetect                      - Initialise the DTMF detection function
 * SDA_DtmfDetect                      - Detect the DTMF signals (uses Goertzel's algorithm)
 * SDA_DtmfDetectAndValidate           - Detect the DTMF signal - as above but validates the output
 * SUF_AsciiToKeyCode                  - Convert ASCII key codes to keypad codes
 * SUF_KeyCodeToAscii                  - Convert keypad codes to ASCII

# Speech Processing Functions


 * SIF_PreEmphasisFilter               - Initialise the pre-emphasis filter for the speech signal
 * SDA_PreEmphasisFilter               - Apply the pre-emphasis filter to the speech signal
 * SIF_DeEmphasisFilter                - Initialise the de-emphasis filter for the speech signal
 * SDA_DeEmphasisFilter                - Apply the de-emphasis filter to the speech signal
 * SDA_AdpcmEncoder                    - Compress the speech using Adaptive Differential Pulse Coded Modulation (ADPCM)
 * SDA_AdpcmEncoderDebug               - Compress the speech using Adaptive Differential Pulse Coded Modulation (ADPCM) - with debug output
 * SDA_AdpcmDecoder                    - Decompress the speech using Adaptive Differential Pulse Coded Modulation (ADPCM)

# Maths / Minimum and Maximum Functions

# Minimum / Maximum Functions
 * SDA_Max                             - Find the maximum value in a array
 * SDA_AbsMax                          - Find the maximum of the absolute values in a array
 * SDA_Min                             - Find the minimum value in a array
 * SDA_AbsMin                          - Find the minimum of the absolute values in a array
 * SDA_Middle                          - Find the middle value in a array
 * SDA_Range                           - Find the range of the minimum and maximum values in a array
 * SDA_MaxPos                          - Find the position of the maximum value in a array
 * SDA_AbsMaxPos                       - Find the position of the maximum of the absolute values in a array
 * SDA_MinPos                          - Find the position of the minimum value in a array
 * SDA_AbsMinPos                       - Find the position of the minimum of the absolute values in a array
 * SDS_Max                             - Return the maximum of 2 numbers
 * SDS_AbsMax                          - Return the maximum of the absolute values of 2 numbers
 * SDS_Min                             - Return the minimum of 2 numbers
 * SDS_AbsMin                          - Return the minimum of 2 the absolute values of numbers
 * SDA_LocalMax                        - Detect the local maximum within an array
 * SDA_LocalAbsMax                     - Detect the local maximum within an array
 * SDA_LocalMin                        - Detect the local maximum within an array
 * SDA_LocalAbsMin                     - Detect the local maximum within an array
 * SDA_Max2                            - Find the maximum value in a array and store in a second
 * SDA_AbsMax2                         - Find the maximum of the absolute values in a array and store in a second
 * SDA_SignedAbsMax2                   - Find the maximum of the absolute values in a array and store in a second, including sign
 * SDA_Min2                            - Find the minimum value in a array and store in a second
 * SDA_AbsMin2                         - Find the minimum of the absolute values in a array and store in a second
 * SDA_SignedAbsMin2                   - Find the minimum of the absolute values in a array and store in a second, including sign
 * SDA_PeakHold                        - Calculate a per array peak hold
 * SDA_PeakHoldPerSample               - Calculate a per sample peak hold
 * SDA_DetectFirstPeakOverThreshold    - Detect the first peak value over a threshold value
 * SDS_Round                           - Round the sample to the nearest integer
 * SDA_Round                           - Round the data to the nearest integer
 * SDS_Clip                            - Clip the sample to an absolute maximum
 * SDA_Clip                            - Clip the samples in the array to an absolute maximum
 * SDS_Threshold                       - Threshold the sample
 * SDA_Threshold                       - Threshold the samples in the array
 * SDS_SoftThreshold                   - Soft threshold the sample
 * SDA_SoftThreshold                   - Soft threshold the samples in the array
 * SDS_ThresholdAndClamp               - Threshold the sample and clamp the results
 * SDA_ThresholdAndClamp               - Threshold the samples in the array and clamp the results
 * SDS_TestOverThreshold               - Test if the sample magnitude is over threshold
 * SDA_TestOverThreshold               - Test if the samples in the array magnitude is over threshold
 * SDS_TestAbsOverThreshold            - Test if absolute value of sample magnitude is over threshold
 * SDA_TestAbsOverThreshold            - Test if absolute value of the samples in the array magnitude is over threshold
 * SDS_Clamp                           - Clamp the sample within limits
 * SDA_Clamp                           - Clamp the samples in the array within limits
 * SDA_SelectMax                       - Select the maximum values from two arrays
 * SDA_SelectMin                       - Select the minimum values from two arrays
 * SDA_SelectMagnitudeSquaredMax       - Select the maximum magnitude squared values from two arrays
 * SDA_SelectMagnitudeSquaredMin       - Select the minimum magnitude squared values from two arrays
 * SDS_SetMinValue                     - Set the minimum value in the output
 * SDA_SetMinValue                     - Set the minimum value in the output
 * SDA_PeakToAverageRatio              - Calculate the peak to average ratio of the signal
 * SDA_PeakToAveragePowerRatio         - Calculate the peak to average power ratio of the signal
 * SDA_PeakToAveragePowerRatioDB       - Calculate the peak to average power ratio of the signal in dB
 * SDA_PeakToAverageRatioComplex       - Calculate the peak to average ratio of the complex signal
 * SDA_PeakToAveragePowerRatioComplex  - Calculate the peak to average power ratio of the complex signal
 * SDA_PeakToAveragePowerRatioComplexDB - Calculate the peak to average power ratio of the complex signal in dB
 * SDA_MovePeakTowardsDeadband         - Shift the data in an array towards a deadband
 * SIF_Envelope                        - Initialize the envelope detector functions
 * SDS_Envelope                        - Envelope detector - per sample
 * SDA_Envelope                        - Envelope detector - on samples in the array
 * SIF_EnvelopeRMS                     - Initialize the RMS energy detection envelope detector functions
 * SDS_EnvelopeRMS                     - RMS energy detection envelope detector - per sample
 * SDA_EnvelopeRMS                     - RMS energy detection envelope detector - on array
 * SIF_EnvelopeHilbert                 - Initialize the Hilbert energy detection envelope detector functions
 * SDS_EnvelopeHilbert                 - Hilbert energy detection envelope detector - per sample
 * SDA_EnvelopeHilbert                 - Hilbert energy detection envelope detector - on array
 * SDS_InterpolateThreePointQuadraticVertexMagnitude          - Quadratic interpolate vertex magnitude (positive and negative) from 3 sample points
 * SDS_InterpolateThreePointQuadraticVertexLocation           - Quadratic interpolate vertex location (positive and negative) from 3 sample points
 * SDS_InterpolateArbitraryThreePointQuadraticVertexMagnitude - Quadratic interpolate vertex magnitude (positive and negative) from 3 sample arbitrary points
 * SDS_InterpolateArbitraryThreePointQuadraticVertexLocation  - Quadratic interpolate vertex location (positive and negative) from 3 sample arbitrary points
 * SDA_InterpolateThreePointQuadraticVertexMagnitude          - Quadratic interpolate vertex magnitude (positive and negative) from 3 points in an array
 * SDA_InterpolateThreePointQuadraticVertexLocation           - Quadratic interpolate vertex location (positive and negative) from 3 points in an array
 * SDA_InterpolateArbitraryThreePointQuadraticVertexMagnitude - Quadratic interpolate vertex magnitude (positive and negative) from 3 arbitrary points in an array
 * SDA_InterpolateArbitraryThreePointQuadraticVertexLocation  - Quadratic interpolate vertex location (positive and negative) from 3 arbitrary points in an array
 * SDA_InterpolateArbitraryThreePointQuadraticPeakVertexMagnitude  - Quadratic interpolate peak vertex magnitude (positive and negative) from 3 arbitrary points in an array
 * SDA_InterpolateArbitraryThreePointQuadraticPeakVertexLocation   - Quadratic interpolate peak vertex location (positive and negative) from 3 arbitrary points in an array
 * SDA_FirstMinVertex                  - Find first minimum virtex in the array
 * SDA_FirstMinVertexPos               - Find first maximum virtex in the array
 * SDA_FirstMaxVertex                  - Find first minimum virtex in the array
 * SDA_FirstMaxVertexPos               - Find first maximum virtex in the array
 * SDA_NLargest                        - Save the N largest values and set the remainder to zero
 * SDA_NSmallest                       - Save the N smallest values and set the remainder to zero

 ## Math Functions

 * SDA_Divide                          - Divide a vector array by a scalar
 * SDA_Divide2                         - Divide a vector array by another vector
 * SDA_Multiply                        - Multiply a vector array by a scalar
 * SDA_Multiply2                       - Multiply two vector arrays
 * SDS_ComplexMultiply                 - Multiply one complex vector by another
 * SDS_ComplexInverse                  - Invert a complex vector
 * SDA_ComplexInverse                  - Invert a complex vectors in the array
 * SDS_ComplexDivide                   - Divide one complex vector by another
 * SDA_ComplexScalarMultiply           - Multiply a complex vector array by a scalar
 * SDA_ComplexMultiply2                - Multiply the complex vectors in an array by the complex vectors in another array
 * SDA_ComplexScalarDivide             - Divide a complex vector array by a scalar
 * SDA_ComplexDivide2                  - Complex divide the elements in one array by the elements in another
 * SDA_RealDotProduct                  - Dot product of two real vectors
 * SDA_ComplexDotProduct               - Dot product of two complex vectors
 * SDA_ComplexScalarMultiply           - Multiply a complex vector array by a scalar
 * SDA_Add_N                           - Add N arrays of data
 * SDA_WeightedSum                     - Weighted sum of two arrays of data
 * SDA_Subtract2                       - Subtract the contents of one array from another
 * SDA_Offset                          - Offset a signal by a scalar amount
 * SDA_PositiveOffset                  - Offset the data to ensure that all the values are positive and the smallest value is zero
 * SDA_NegativeOffset                  - Offset the data to ensure that all the values are negative and the largest value is zero
 * SDA_Negate                          - Negate all the entries in a array
 * SDA_Inverse                         - Return reciprocal of all the values in the array
 * SDA_Square                          - Return the square of all the values in the array
 * SDA_Sqrt                            - Return the square root of all the values in the array
 * SDA_Difference                      - Returns the difference (always positive) between the data in the two arrays
 * SDA_SumOfDifferences                - Returns the sum of the differences (always positive) between the data in the two arrays
 * SDS_Roots                           - Returns the real roots of the equation: ax^2 + bx + c = 0
 * SDS_Factorial                       - Returns the factorial of the real input value
 * SDS_Permutations                    - Returns the number of permutations for a data set
 * SDS_Combinations                    - Returns the number of combinations for a data set
 * SIF_OverlapAndAddLinear             - Initialise the linear overlap and add functions
 * SDA_OverlapAndAddLinear             - Overlap the data in two arrays and add the overlapped data values
 * SDA_OverlapAndAddLinearWithClip     - Overlap the data in two arrays and add the overlapped data values - Clip the sum to avoid overflow
 * SDA_OverlapAndAddArbitrary          - Overlap the data in two arrays and add the overlapped data values - Adds an arbitrary scaling function to the overlapped data
 * SDA_OverlapAndAddArbitraryWithClip  - Overlap the data in two arrays and add the overlapped data values -
 * Adds an arbitrary scaling function to the overlapped data and clip the sum to avoid overflow
 * SDS_DegreesToRadians                - Convert a value from degrees to radians
 * SDA_DegreesToRadians                - Convert an array of values from degrees to radians
 * SDS_RadiansToDegrees                - Convert a value from radians to degrees
 * SDA_RadiansToDegrees                - Convert an array of values from radians to degrees
 * SDS_DetectNAN                       - Test if the sample is NAN or INF
 * SDA_DetectNAN                       - Test if any of the samples in the array are NAN or INF


# DSP utility Functions


 * SDA_Rotate                          - Rotate the data in a array by N samples
 * SDA_Reverse                         - Reverse the data in a array
 * SDA_Scale                           - Scale an arrays contents to a maximum value
 * SDA_MeanSquare                      - Calculate the MS value of the data
 * SDA_MeanSquareError                 - Calculate the MS error value of the data
 * SDA_RootMeanSquare                  - Calculate RMS value of the data
 * SDA_Magnitude                       - Calculate RMS magnitude from complex arrays
 * SDA_MagnitudeSquared                - Calculate mean squared magnitude from complex arrays
 * SDS_Magnitude                       - Calculate RMS magnitude from complex sample
 * SDS_MagnitudeSquared                - Calculate mean squared magnitude from a single complex input
 * SDS_Phase                           - Calculate phase from a complex sample
 * SDA_Phase                           - Calculate phase from a complex array
 * SDA_PhaseWrapped                    - Return the wrapped phase of the signal around +/- PI
 * SDA_PhaseUnwrapped                  - Return the unwrapped of the signal
 * SDA_MagnitudeAndPhaseWrapped        - Return the magnitude and the phase wrapped of the signal
 * SDA_MagnitudeAndPhaseUnWrapped      - Return the magnitude and the phase unwrapped of the signal
 * SDA_MagnitudeSquaredAndPhaseWrapped - Return the magnitude and the phase wrapped of the signal
 * SDA_MagnitudeSquaredAndPhaseUnWrapped - Return the magnitude and the phase unwrapped of the signal
 * SDA_PhaseWrap                       - Wrap the phase of the signal around +/- PI
 * SDA_PhaseUnwrap                     - Unwrap the phase of the signal
 * SDS_Log2                            - Return the logartitm to base 2 of the sample
 * SDA_Log2                            - Return the logartitm to base 2 for all the samples in the array
 * SDS_LogN                            - Return the logartitm to base N of the sample
 * SDA_LogN                            - Return the logartitm to base N for all the samples in the array
 * SDA_Copy                            - Duplicate the contents of one array in another
 * SDA_CopyWithStride                  - Duplicate the contents with variable indexing
 * SIF_CopyWithOverlap                 - Initilize the function for copying two arrays to a third with some overlap of the input arrays
 * SDA_CopyWithOverlap                 - Copy two arrays to a third with some overlap of the input arrays
 * SIF_CopyWithIndex                   - Initialize copy with index function
 * SDA_CopyWithIndex                   - Copy contents of one array into a second, with increasing index into the first array, on each copy
 * SDA_RectangularToPolar              - Convert rectangular data to polar
 * SDA_PolarToRectangular              - Convert polar data to rectangular
 * SDA_20log10                         - Calculate 20 * log base 10 of vector values on an array basis
 * SDA_10log10                         - Calculate 10 * log base 10 of vector values on an array basis
 * SDS_20log10                         - Calculate 20 * log base 10 of vector values on a per sample basis
 * SDS_10log10                         - Calculate 10 * log base 10 of vector values on a per sample basis
 * SDA_LogMagnitude                    - Calculate the log magnitude of the vector
 * SDA_LogMagnitudeAndPhaseWrapped     - Return the log magnitude and the phase wrapped of the signal
 * SDA_LogMagnitudeAndPhaseUnWrapped   - Return the log magnitude and the phase unwrapped of the signal
 * SDA_Lengthen                        - Zero pad a vector to another length
 * SDA_Shorten                         - Truncate a vector to a new length
 * SIF_ReSize                          - Initialize the array resize function
 * SDA_ReSize                          - Resize the array
 * SDA_ReSizeInput                     - Resize the input
 * SDA_ReSizeOutput                    - Resize the output
 * SDA_Fill                            - Fill a array with a scalar value
 * SDA_Clear                           - Clear the contents of a array to 0.0
 * SDA_Histogram                       - Produce a histogram of an arrays data
 * SDA_HistogramCumulative             - Produce a cumulative histogram of an arrays data
 * SDA_HistogramExtended               - Produce an extended histogram of an arrays data
 * SDA_HistogramExtendedCumulative     - Produce an extended cumulative histogram of an arrays data
 * SIF_Histogram                       - Initilize the histogram functions
 * SDA_HistogramEqualize               - Equalize the histogram of the data
 * SDA_Quantize                        - Quantize the data in the array to N bits
 * SDS_Quantize                        - Quantize the saample to N bits
 * SDA_Quantize_N                      - Quantize, to integer N, the data in the array to N bits
 * SDS_Quantize_N                      - Quantize, to integer N, the saample to N bits
 * SDA_Abs                             - Calculate the absolute values in an array
 * SDS_PeakValueToBits                 - Convert the number of bits in an integer to the peak value possible
 * SDS_BitsToPeakValue                 - For a given value, provide the number of bits required to represent it
 * SDS_LinearTodBm                     - Convert a linear sample to dBm
 * SDA_LinearTodBm                     - Convert a linear array of samples to dBm
 * SDS_dBmToLinear                     - Convert a dBm sample to linear
 * SDA_dBmToLinear                     - Convert a dBm array of samples to linear
 * SDS_Compare                         - Compare the contents of two samples
 * SDA_Compare                         - Compare the contents of two arrays
 * SDS_CompareComplex                  - Compare the contents of two complex samples
 * SDA_CompareComplex                  - Compare the contents of two arrays of complex numbers
 * SDS_Int                             - Integer value of the sample
 * SDS_Frac                            - Fractional value of the sample
 * SDS_AbsFrac                         - Absolute value of the fractional value of the sample
 * SDA_Int                             - Integer value of the array samples
 * SDA_Frac                            - Fractional value of the array samples
 * SDA_AbsFrac                         - Absolute value of the fractional value of the array samples
 * SDA_SetRange                        - Set the range of the array samples
 * SDA_SetMean                         - Set the mean of the array samples
 * SDA_RealSpectralInverse             - Spectrum inv. on real time domain data
 * SDA_ComplexSpectralInverse          - Spectrum inv. on Complex time domain data
 * SDA_FdInterpolate                   - Interpolate a spectrum, to change pitch
 * SDA_FdInterpolate2                  - Interpolate a spectrum, to change pitch
 * SDS_TdPitchShift                    - Pitch shift a signal in the time domain - per sample
 * SDA_TdPitchShift                    - Pitch shift a signal in the time domain - array
 * SDS_EchoGenerate                    - Superimpose echo and reverb on a signal
 * SDA_Power                           - Raise all the entries in a array to a power
 * SDS_Polynomial                      - Evaluate the polynomial on the sample
 * SDA_Polynomial                      - Evaluate the polynomial on the data
 * SDS_Modulo                          - Rewrite the data as modulo N data
 * SDA_Modulo                          - Rewrite the array of data as modulo N data
 * SDA_AgcPeak                         - Control the gain using the peak level
 * SIF_AgcMeanAbs                      - Initialize the AGC to control the gain using the mean absolute value
 * SDA_AgcMeanAbs                      - Control the gain using the mean absolute value
 * SIF_AgcMeanSquared                  - Initialize the AGC to control the gain using the mean squared value
 * SDA_AgcMeanSquared                  - Control the gain using the mean squared value
 * SDA_GroupDelay                      - Return the group delay of the phase signal
 * SDA_ZeroCrossingDetect              - Find zero crossings in an array
 * SDS_ZeroCrossingDetect              - Find zero crossings on a per sample basis
 * SDA_FirstZeroCrossingLocation       - Find the location of the first zero crossing on an array basis
 * SDA_ZeroCrossingCount               - Count the number of zero crossings on an array basis
 * SDA_LevelCrossingDetect             - Find level crossings in an array
 * SDS_LevelCrossingDetect             - Find level crossings on a per sample basis
 * SDA_FirstLevelCrossingLocation      - Find the location of the first level crossing on an array basis
 * SDA_LevelCrossingCount              - Count the number of level crossings on an array basis
 * SDA_Trigger                         - Oscilloscope style trigger function
 * SDA_ClearLocation                   - Set the value at a location to zero
 * SDA_SetLocation                     - Set the value at a location to the given value
 * SDA_SortMinToMax                    - Sort the values in a array - minimum first
 * SDA_SortMaxToMin                    - Sort the values in a array - maximum first
 * SDA_SortMinToMax2                   - Sort the values in two arrays - minimum first
 * SDA_SortMaxToMin2                   - Sort the values in two arrays - maximum first
 * SDA_SortIndexed                     - Sort the values in a array - indexed
 * SDS_CountOneBits                    - Count the number of '1' bits in the data value
 * SDS_CountZeroBits                   - Count the number of '0' bits in the data value
 * SDS_CountLeadingOneBits             - Count the number of leading '1' bits in the data value
 * SDS_CountLeadingZeroBits            - Count the number of leading '0' bits in the data value
 * SDA_Sign                            - Return the sign of all the samples in the array
 * SDA_Swap                            - Swap the order of the data in the array
 * SUF_ModuloIncrement                 - Increment all of the values in the source array using modulo arithmatic
 * SUF_ModuloDecrement                 - Decrement all of the values in the source array using modulo arithmatic
 * SUF_IndexModuloIncrement            - Increment all of the index values in the source array using modulo arithmatic
 * SUF_IndexModuloDecrement            - Decrement all of the index values in the source array using modulo arithmatic
 * SDA_Find                            - Locate all the values in the source array that match the specification
 * SDA_FindValue                       - Locate all the values in the source array that match the value
 * SIF_DeGlitch                        - Initialize the de-glitch functions
 * SDS_DeGlitch                        - Deglitch the data samples
 * SDA_DeGlitch                        - Deglitch the samples in the array
 * SDA_RemoveDuplicates                - Remove duplicate values in the array
 * SDA_FindAllDuplicates               - Find all duplicate values in the array
 * SDA_FindFirstDuplicates             - Find first duplicates in the array
 * SDA_FindSortAllDuplicates           - Find and sort, smallest to largest, all duplicate values in the array
 * SDA_FindSortFirstDuplicates         - Find and sort, smallest to largest, first duplicates in the array
 * SDA_Shuffle                         - Shuffle the data in the array
 * SDA_SigLibDataToFix                 - Convert the data from native SigLib to fixed point format
 * SDA_FixToSigLibData                 - Convert the data from fixed point format to native SigLib
 * SDA_SigLibDataToImageData           - Convert the data from native SigLib to SigLib fixed point image format
 * SDA_ImageDataToSigLibData           - Convert the data from SigLib fixed point image to native SigLib format
 * SDA_SigLibDataToFix16               - Convert the data from native SigLib to 16 bit fixed point format
 * SDA_Fix16ToSigLibData               - Convert the data from 16 bit fixed point format to native SigLib
 * SDA_SigLibDataToFix32               - Convert the data from native SigLib to 32 bit fixed point format
 * SDA_Fix32ToSigLibData               - Convert the data from 32 bit fixed point format to native SigLib
 * SDS_SigLibDataToQFormatInteger      - Convert the data from native SigLib to Q format - per sample basis
 * SDS_QFormatIntegerToSigLibData      - Convert the data from Q format to native SigLib data format - per sample basis
 * SDA_SigLibDataToQFormatInteger      - Convert the data from native SigLib to Q format - array basis
 * SDA_QFormatIntegerToSigLibData      - Convert the data from Q format to native SigLib data format - array basis


# Control Functions
 * SDS_Pid                             - Apply a pid control loop
 * SDA_Pwm                             - Pulse Width Modulation function


# Order Analysis Functions

 * SDA_ExtractOrder                    - Extract the orders in a data set
 * SDA_SumLevel                        - Sum the levels in a data set
 * SDA_SumLevelWholeSpectrum           - Sum the levels over a whole spectrum
 * SIF_OrderAnalysis                   - Initialize the order analysis functions
 * SDA_OrderAnalysis                   - Perform order analysis

# Statistical Analysis Functions

 * SDA_Sum                             - Sum all the entries in the array
 * SDA_AbsSum                          - Sum the absolute values of all the entries in an array
 * SDA_SumOfSquares                    - Return the sum of squares of all the entries in the array
 * SDA_Mean                            - Calculate the arithmetic mean (average) of the array
 * SDA_AbsMean                         - Calculate the absolute value of the arithmetic mean (average) of the array
 * SDA_SubtractMean                    - Subtract the arithmetic mean (average) of the values in the array from all the values
 * SDA_SubtractMax                     - Subtract the maximum value in the array from all the values
 * SDA_SampleSd                        - Calculate the sample standard deviation of the array
 * SDA_PopulationSd                    - Calculate the population standard deviation of the array
 * SDA_UnbiasedVariance                - Calculate the unbiased variance of the array
 * SDA_Median                          - Calculate the median value of the array
 * Top

# Regression Analysis Functions

 * SDA_LinraConstantCoeff              - Return lin. regression const. coeff.
 * SDA_LinraRegressionCoeff            - Return lin. regression regress. coeff.
 * SDA_LinraCorrelationCoeff           - Return lin. regression correlat. coeff.
 * SDA_LinraEstimateX                  - Return lin. regression X estimate
 * SDA_LinraEstimateY                  - Return lin. regression Y estimate
 * SDA_LograConstantCoeff              - Return log. regression const. coeff.
 * SDA_LograRegressionCoeff            - Return log. regression regress. coeff.
 * SDA_LograCorrelationCoeff           - Return log. regression correlat. coeff.
 * SDA_LograEstimateX                  - Return log. regression X estimate
 * SDA_LograEstimateY                  - Return log. regression Y estimate
 * SDA_ExpraConstantCoeff              - Return expon. regression const. coeff.
 * SDA_ExpraRegressionCoeff            - Return expon. regression regress. coeff.
 * SDA_ExpraCorrelationCoeff           - Return expon. regression correlat. coeff.
 * SDA_ExpraEstimateX                  - Return expon. regression X estimate
 * SDA_ExpraEstimateY                  - Return expon. regression Y estimate
 * SDA_PowraConstantCoeff              - Return power regression const. coeff.
 * SDA_PowraRegressionCoeff            - Return power regression regress. coeff.
 * SDA_PowraCorrelationCoeff           - Return power regression correlat. coeff.
 * SDA_PowraEstimateX                  - Return power regression X estimate
 * SDA_PowraEstimateY                  - Return power regression Y estimate
 * SDA_Detrend                         - Detrend the data in the array
 * SDA_ExtractTrend                    - Extract the trend in the data

# Trigonometric Functions

 * SDA_Sin                             - Return the sines of the data
 * SDA_Cos                             - Return the cosines of the data
 * SDA_Tan                             - Return the tangents of the data
 * SIF_FastSin                         - Initialise fast sine transform (look up table) function
 * SDA_FastSin                         - Perform fast sine (look up table) function on an array basis
 * SDS_FastSin                         - Perform fast sine (look up table) function on a per sample basis
 * SIF_FastCos                         - Initialise fast cosine transform (look up table) function
 * SDA_FastCos                         - Perform fast cosine (look up table) function on an array basis
 * SDS_FastCos                         - Perform fast cosine (look up table) function on a per sample basis
 * SIF_FastSinCos                      - Initialise fast sine and cosine transform (look up table) function
 * SDA_FastSinCos                      - Perform fast sine and cosine (look up table) function on an array basis
 * SDS_FastSinCos                      - Perform fast sine and cosine (look up table) function on a per sample basis
 * SIF_FastTan                         - Initialise fast tangent transform (look up table) function
 * SDA_FastTan                         - Perform fast tangent (look up table) function on an array basis
 * SDS_FastTan                         - Perform fast tangent (look up table) function on a per sample basis
 * SIF_QuickSin                        - Initialise quick sine transform (look up table) function
 * SDA_QuickSin                        - Perform quick sine (look up table) function on an array basis
 * SDS_QuickSin                        - Perform quick sine (look up table) function on a per sample basis
 * SIF_QuickCos                        - Initialise quick cosine transform (look up table) function
 * SDA_QuickCos                        - Perform quick cosine (look up table) function on an array basis
 * SDS_QuickCos                        - Perform quick cosine (look up table) function on a per sample basis
 * SIF_QuickSinCos                     - Initialise quick sine and cosine transform (look up table) function
 * SDA_QuickSinCos                     - Perform quick sine and cosine (look up table) function on an array basis
 * SDS_QuickSinCos                     - Perform quick sine and cosine (look up table) function on a per sample basis
 * SIF_QuickTan                        - Initialise quick tangent transform (look up table) function
 * SDA_QuickTan                        - Perform quick tangent (look up table) function on an array basis
 * SDS_QuickTan                        - Perform quick tangent (look up table) function on a per sample basis
 * SDA_Sinc                            - Calculate sin(x)/x sinc function on an array
 * SDS_Sinc                            - Calculate sin(x)/x sinc function on a per sample basis
 * SIF_QuickSinc                       - Initialize the function to calculate the quick sin(x)/x sinc (look up table) function
 * SDA_QuickSinc                       - Calculate the quick sin(x)/x sinc (look up table) function on an array
 * SDS_QuickSinc                       - Calculate the quick sin(x)/x sinc (look up table) function on a per sample basis

# Complex Vector Functions


 * SCV_Real                            - Extract the real part of a complex vector
 * SCV_Imaginary                       - Extract the imaginary part of a complex vector
 * SCV_Polar                           - Create a complex polar vector from real data
 * SCV_Rectangular                     - Create a complex rectangular vector from real data
 * SCV_PolarToRectangular              - Convert polar data to rectangular
 * SCV_RectangularToPolar              - Convert rectangular data to polar
 * SCV_Sqrt                            - Square root a complex vector
 * SCV_Inverse                         - Invert a complex vector
 * SCV_Conjugate                       - Return the complex conjugate of a vector
 * SCV_Magnitude                       - Return the real absolute magnitude of vector
 * SCV_MagnitudeSquared                - Return the square of the real absolute magnitude of vector
 * SCV_Phase                           - Return the phase of vector
 * SCV_Multiply                        - Multiply two complex vectors
 * SCV_Divide                          - Divide one complex vector by another
 * SCV_Add                             - Add two complex vectors
 * SCV_Subtract                        - Return the difference between two vectors
 * SCV_Log                             - Return log of a complex vector
 * SCV_Exp                             - Return exponent of a complex vector
 * SCV_Expj                            - Return the complex exponential of the real input (e^jTheta)
 * SCV_Pow                             - Raise complex vector to a real power
 * SCV_VectorAddScalar                 - Add a scalar to a vector
 * SCV_VectorSubtractScalar            - Subtract a scalar from a vector
 * SCV_VectorMultiplyScalar            - Multiply a vector by a scalar
 * SCV_VectorDivideScalar              - Divide a vector by a scalar
 * SCV_ScalarSubtractVector            - Subtract a vector from a scalar
 * SCV_Roots                           - Returns the complex roots of the equation: ax^2 + bx + c = 0
 * SCV_Copy                            - Copy one vector to another
 * SCV_Compare                         - Compare two vectors
 * SDA_CreateComplexRect               - Create a rectangular complex array (interleaved real and imaginary values)
 * SDA_CreateComplexPolar              - Create a polar complex array (interleaved magnitude and phase values)
 * SDA_ExtractComplexRect              - Separate a rectangular complex array (interleaved real and imaginary values)
 * SDA_ExtractComplexPolar             - Separate a polar complex array (interleaved magnitude and phase values)
 * SDA_ClearComplexRect                - Clear the contents of a rectangular complex array to 0+j0
 * SDA_ClearComplexPolar               - Clear the contents of a polar complex array to 0+j0
 * SDA_FillComplexRect                 - Fill the contents of a rectangular complex array to a constant complex value
 * SDA_FillComplexPolar                - Fill the contents of a polar complex array to a constant complex value
 * SDA_ComplexRectangularToPolar       - Convert complex rectangular numbers to polar
 * SDA_ComplexPolarToRectangular       - Convert complex polar numbers to rectangular
 * SDA_RectangularToPolar              - Convert rectangular numbers to polar
 * SDA_PolarToRectangular              - Convert polar numbers to rectangular
 * SDA_ComplexRectSqrt                 - Square root of all the values in an array of complex numbers
 * SDA_ComplexRectInverse              - Inverse of all the values in an array of complex numbers
 * SDA_ComplexRectConjugate            - Complex conjugate of all the values in an array of complex numbers
 * SDA_ComplexRectMagnitude            - Magnitude of all the values in an array of complex numbers
 * SDA_ComplexRectMagnitudeSquared     - Magnitude squared of all the values in an array of complex numbers
 * SDA_ComplexRectPhase                - Phase of all the values in an array of complex numbers
 * SDA_ComplexRectMultiply             - Multiply all the values in two arrays of complex numbers
 * SDA_ComplexRectDivide               - Divide all the values in two arrays of complex numbers
 * SDA_ComplexRectAdd                  - Add all the values in two arrays of complex numbers
 * SDA_ComplexRectSubtract             - Subtract all the values in two arrays of complex numbers
 * SDA_ComplexRectLog                  - Log of all the values in an array of complex numbers
 * SDA_ComplexRectExp                  - Exponential of all the values in an array of complex numbers
 * SDA_ComplexRectExpj                 - Complex exponential of all the values in an array of complex numbers
 * SDA_ComplexRectPow                  - Raise of all the values in an array of complex numbers to a power
 * SDA_ComplexRectVectorAddScalar      - Add a scalar value to all the values in an array of complex numbers
 * SDA_ComplexRectVectorSubtractScalar - Subtract a scalar value from all the values in an array of complex numbers
 * SDA_ComplexRectVectorMultiplyScalar - Multiply all the values in an array of complex numbers by a scalar
 * SDA_ComplexRectVectorDivideScalar   - Divide all the values in an array of complex numbers by a scalar
 * SDA_ComplexRectScalarSubtractVector - Subtract all the values in an array of complex numbers from a scalar number
 * SDA_ComplexRectLinearInterpolate    - Interpolate the values in an array of complex numbers in a recgtangular mode
 * SDA_ComplexPolaInterpolate          - Interpolate the values in an array of complex numbers in a polar mode


# Two Dimensional Matrix Functions


 * SMX_Transpose                       - Transpose a two dimensional matrix - also know as 'corner turn'
 * SMX_Multiply                        - Multiply two two dimensional matrices
 * SMX_CreateIdentity                  - Create an identity matrix
 * SMX_Inverse2x2                      - Invert a 2x2 matrix
 * SMX_ComplexInverse2x2               - Invert a 2x2 matrix of complex numbers
 * SMX_Inverse                         - Invert a matrix, using Crout's LU reduction algorithm
 * SMX_LuDecompose                     - LU decompose a matrix
 * SMX_LuSolve                         - Solve for an LU decomposed matrix
 * SMX_Determinant                     - Calculate the determinant of a matrix
 * SMX_LuDeterminant                   - Calculate the determinant of an LU decomposed matrix
 * SMX_RotateClockwise                 - Rotate all of the matrix values in a clockwise direction
 * SMX_RotateAntiClockwise             - Rotate all of the matrix values in an anti-clockwise direction
 * SMX_Reflect                         - Reflect the values in a matrix about a vertical axis
 * SMX_Flip                            - Reflect the values in a matrix about a horizontal axis
 * SMX_InsertRow                       - Insert a row in a matrix
 * SMX_ExtractRow                      - Extract a row from a matrix
 * SMX_InsertColumn                    - Insert a column in a matrix
 * SMX_ExtractColumn                   - Extract a column from a matrix
 * SMX_InsertNewRow                    - Insert a new row into a matrix - shift the rows below the new one down one level
 * SMX_DeleteOldRow                    - Extract a row from a matrix - shift the rows below the extracted one up one level
 * SMX_InsertNewColumn                 - Insert a new column into a matrix - shift the columns to the right of the new one to the right one level
 * SMX_DeleteOldColumn                 - Extract a column from a matrix - shift the columns to the right of the extracted one to the left one level
 * SMX_InsertRegion                    - Insert a region into a matrix
 * SMX_ExtractRegion                   - Extract a region from a matrix
 * SMX_InsertDiagonal                  - Insert a diagonal into a matrix
 * SMX_ExtractDiagonal                 - Extract a diagonal from a matrix
 * SMX_SwapRows                        - Swap two rows in a matrix
 * SMX_SwapColumns                     - Swap two columns in a matrix
 * SMX_Sum                             - Add two matrices together
 * SMX_ShuffleColumns                  - Shuffle columns in a matrix
 * SMX_ShuffleRows                     - Shuffle rows in a matrix
 * SMX_ExtractCategoricalColumn        - Extract categorical (last) column from a matrix
 * SMX_Copy                            - Copy a two dimensional matrix
 * SMX_Add                             - Add two two dimensional matrices
 * SMX_Subtract                        - Subtract one two D. matrix from another
 * SMX_MultiplyPiecewise               - Piecewise multiply two two D. matrices
 * SMX_ScalarMultiply                  - Multiply a two D. matrix by a scalar

# File I/O Functions


 * SUF_BinReadData                     - Read data from a binary file
 * SUF_BinWriteData                    - Write data to a binary file
 * SUF_BinReadFile                     - Read a complete binary file
 * SUF_BinWriteFile                    - Write a complete binary file
 * SUF_CsvReadData                     - Read a complete csv file
 * SUF_CsvWriteFile                    - Write a complete csv file
 * SUF_CsvReadFile                     - Read data from a csv file
 * SUF_CsvWriteData                    - Write data to a csv file
 * SUF_CsvReadMatrix                   - Read matrix from a csv file
 * SUF_DatReadData                     - Read data from a dat file
 * SUF_DatWriteData                    - Write data to a dat file
 * SUF_DatReadHeader                   - Read header from a dat file
 * SUF_DatWriteHeader                  - Write header to a binary file
 * SUF_SigReadData                     - Read data from a sig file
 * SUF_SigWriteData                    - Write data to a sig file
 * SUF_SigReadFile                     - Read a complete sig file
 * SUF_SigWriteFile                    - Write a complete sig file
 * SUF_WavReadData                     - Read data from a wav file
 * SUF_WavWriteData                    - Write data to a wav file
 * SUF_WavReadWord                     - Read word from a wav file
 * SUF_WavReadInt                      - Read integer from a wav file
 * SUF_WavWriteWord                    - Write word to a wav file
 * SUF_WavWriteInt                     - Write integer to a wav file
 * SUF_WavReadHeader                   - Read header from a wav file
 * SUF_WavWriteHeader                  - Write header to a wav file
 * SUF_WavDisplayInfo                  - Display header information from a wav file
 * SUF_WavSetInfo                      - Set header information for a wav file
 * SUF_WavFileLength                   - Write file length to a wav file
 * SUF_WavReadFile                     - Read a complete wav file
 * SUF_WavWriteFile                    - Write a complete wav file
 * SUF_WavWriteFileScaled              - Write scaled data to a complete wav file
 * SUF_XmtReadData                     - Read data from an xmt file
 * SUF_WriteWeightsIntegerCFile        - Write neural network weights to an integer C header file
 * SUF_WriteWeightsFloatCFile          - Write neural network weights to a floating point C header file
 * SUF_WriteWeightsBinaryFile          - Write neural network weights to a binary file
 * SUF_ReadWeightsBinaryFile           - Read neural network weights from a binary file


# Machine Learning Functions

## Machine Learning Functions

 * SDA_TwoLayer2CategoryNetworkFit         - Train 2 layer, 2 category Neural Network
 * SDA_TwoLayer2CategoryNetworkPredict     - Predict category using 2 layer, 2 category Neural Network
 * SDA_TwoLayerNCategoryNetworkFit         - Train 2 layer, arbitrary number of categories Neural Network
 * SDA_TwoLayerNCategoryNetworkPredict     - Predict category using 2 layer, arbitrary number of categories Neural Network
 * SDS_ActivationReLU                      - ReLU activation function on a single sample
 * SDA_ActivationReLU                      - ReLU activation function on an array of samples
 * SDS_ActivationReLUDerivative            - Derivative of ReLU activation function on a single sample
 * SDA_ActivationReLUDerivative            - Derivative of ReLU activation function on an array of samples
 * SDS_ActivationLeakyReLU                 - Leaky ReLU activation function on a single sample
 * SDA_ActivationLeakyReLU                 - Leaky ReLU activation function on an array of samples
 * SDS_ActivationLeakyReLUDerivative       - Derivative of leaky ReLU activation function on a single sample
 * SDA_ActivationLeakyReLUDerivative       - Derivative of leaky ReLU activation function on an array of samples
 * SDS_ActivationLogistic                  - Logistic activation function on a single sample
 * SDA_ActivationLogistic                  - Logistic activation function on an array of samples
 * SDS_ActivationLogisticDerivative        - Derivative of logistic activation function on a single sample
 * SDA_ActivationLogisticDerivative        - Derivative of logistic activation function on an array of samples
 * SDS_ActivationTanH                      - Tan H activation function on a single sample
 * SDA_ActivationTanH                      - Tan H activation function on an array of samples
 * SDS_ActivationTanHDerivative            - Derivative of tan H activation function on a single sample
 * SDA_ActivationTanHDerivative            - Derivative of tan H activation function on an array of samples

# Utility Macros And Functions

## Utility Functions

 * SUF_SiglibVersion                   - Returns the current SigLib version number
 * SUF_Debugvfprintf                   - Varaible parameter version of SUF_Debugfprintf
 * SUF_DebugPrintArray                 - Print the contents of an array to the log file
 * SUF_DebugPrintComplexArray          - Print the contents of a complex array to the log file
 * SUF_DebugPrintMatrix                - Print the contents of a matrix to the log file
 * SUF_DebugPrintPolar                 - Print the value of a polar complex number to the log file
 * SUF_DebugPrintRectangular           - Print the value of a rectangular complex number to the log file
 * SUF_DebugPrintIIRCoefficients       - Print the IIR filter coefficients array to the log file
 * SUF_ClearDebugfprintf               - Clear the debug log file
 * SUF_Debugfprintf                    - Debug fprintf function - prints to a log file
 * SUF_Debugvprintf                    - Debug vfprintf function - prints to a log file
 * SUF_PrintArray                      - Print the contents of an array to the console
 * SUF_PrintComplexArray               - Print the contents of a complex array to the console
 * SUF_PrintMatrix                     - Print the contents of a matrix to the console
 * SUF_PrintPolar                      - Print the value of a polar complex number to the console
 * SUF_PrintRectangular                - Print the value of a rectangular complex number to the console
 * SUF_PrintIIRCoefficients            - Print the IIR filter coefficients array to the console
 * SUF_PrintCount                      - Print a count to keep track of iterations
 * SUF_DebugPrintCount                 - Print a count to keep track of iterations to the log file
 * SUF_PrintHigher                     - Print the higher of two numbers
 * SUF_PrintLower                      - Print the lower of two numbers
 * SUF_DebugPrintHigher                - Print the higher of two numbers to the log file
 * SUF_DebugPrintLower                 - Print the lower of two numbers to the log file
 * SUF_DebugPrintInfo                  - Print SigLib library version information to the log file
 * SUF_DebugPrintLine                  - Print source code line information to the log file
 * SUF_DebugPrintTime                  - Print time information to the log file
 * SUF_StrError                        - Print the SigLib error number
 * SUF_MSDelay                         - Delay the code by N ms

# C/C++ Macros

 * SDS_RoundDown                       - Round the scalar value down to the nearest integer value
 * SDS_RoundUp                         - Round the scalar value up to the nearest integer value
 * SDS_RoundToNearest                  - Round the scalar value to the nearest integer value
 * SDS_Odd                             - Returns true if sample is odd
 * SDS_Even                            - Returns true if sample is even
 * SDS_PowerOfTwo                      - Returns true if size is a power of 2
 * SDS_Absolute                        - Return the absolute value of a number
 * SDS_Sign                            - Return the sign of a number
 * SDS_BitTest                         - Returns 1 if all bits in mask set, 0 otherwise
 * SUF_NumberOfElements                - Returns number of samples in an array
 * SDS_BitTest                         - Test bit pattern
 * SDS_BitMask                         - Bit mask
 * SDS_Swap                            - Swap the contents of 2 variables
 * SDS_SortN                           - Sort the contents of N variables
 * SDA_Operate                         - Perform a standard math operation (+ - * / ) on vectors
 * SAI_RoundDown                       - Round the integer scalar value down to the nearest integer value
 * SAI_RoundUp                         - Round the integer scalar value up to the nearest integer value
 * SAI_RoundToNearest                  - Round the integer scalar value to the nearest integer value
 * SAI_Odd                             - Return true if the integer value is odd
 * SAI_Even                            - Return true if the integer value is even
 * SAI_PowerOfTwo                      - Return the integer value squared
 * SAI_Absolute                        - Return the absolute value of the integer
 * SAI_Sign                            - Return the sign of the integer value
 * SAI_Log2                            - Return the log base 2 of the integer value
 * SAI_Log4                            - Return the log base 4 of the integer value
 * SUF_BinNumberToFrequency            - Convert the FFT bin number to a frequency
 * SUF_BinNumberToFrequency2           - Convert the FFT bin number to a frequency - version 2
 * SUF_FrequencyToBinNumber            - Convert a frequency to an FFT bin number
 * SUF_FrequencyToBinNumber2           - Convert a frequency to an FFT bin number - version 2

 # Memory Allocation Functions

 * SUF_VectorArrayAllocate             - Allocate vector array
 * SUF_FftCoefficientAllocate          - Allocate FFT twiddle factor coefficient array
 * SUF_IirState_arrayAllocate          - Allocate IIR filter state array
 * SUF_IirCoefficientAllocate          - Allocate IIR filter coefficients  array
 * SUF_QamCarrierArrayAllocate         - Allocate QAM carrier array
 * SUF_QpskCarrierArrayAllocate        - Allocate QPSK carrier array
 * SUF_IndexArrayAllocate              - Allocate index array
 * SUF_ComplexRectArrayAllocate        - Allocate a rectangular complex vector array
 * SUF_ComplexPolarArrayAllocate       - Allocate a polar complex vector array



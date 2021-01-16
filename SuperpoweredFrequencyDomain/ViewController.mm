#import "ViewController.h"
#import "SuperpoweredIOSAudioIO.h"
#import "Freq_Transforms.h"
#include "SuperpoweredFrequencyDomain.h"
#include "SuperpoweredSimple.h"
#include "Data.h"

int nFFT=512;
//int FFT=1024;
//static float *magnitudeRight =(float*)calloc(nFFT/2, sizeof(float));
//static float *magnitudeLeft =(float*)calloc(nFFT/2, sizeof(float));
static float *noise_mean = (float*)calloc(nFFT, sizeof(float));
static float *noise_mu2 = (float*)calloc(nFFT, sizeof(float));
static float *sig2 = (float*)calloc(nFFT, sizeof(float));
static float *gammak = (float*)calloc(nFFT, sizeof(float));
static float *ksi = (float*)calloc(nFFT, sizeof(float));
static float *Xk_prev = (float*)calloc(nFFT, sizeof(float));
static float *log_sigma_k = (float*)calloc(nFFT, sizeof(float));
static float *noise_pow = (float*)calloc(nFFT, sizeof(float));
static float *vk = (float*)calloc(nFFT, sizeof(float));
static float *hw = (float*)calloc(nFFT, sizeof(float));
static float *evk = (float*)calloc(nFFT, sizeof(float));
static float *Lambda = (float*)calloc(nFFT, sizeof(float));
static float *pSAP = (float*)calloc(nFFT, sizeof(float));
static float* ensig = ( float *)calloc(nFFT, sizeof( float));
static float* HPF = ( float*)calloc(2*(nFFT), sizeof( float));
static float* H = ( float*)calloc((nFFT), sizeof( float));
static float *fftmagoutput = (float*)calloc(nFFT, sizeof(float));
//long double* nsig = (long double *)calloc(nFFT, sizeof(long double));
float PRT;
float N;
float muu=0.5;
float v=.01;
float aa = 0.98;
float eta = 0.15;
float beta=0.5;
float max;
float ksi_min = (float)pow(10,((float)-25 /(float)10));
float sum_log_sigma_k = 0;
float vad_decision;
float qk = 0.3;
float qkr = (1 - qk) / qk;;
//float epsilon =  (float)pow(8.854,-12);
float epsilon = 0.001;
float total_noisepower=0;
float total_speechpower=0;
//float SNR_db;
char SNR_db_char;
int count=0;
int SPU = 0;
static float SNR_db;
//float *magnitudeLeft, *magnitudeRight, *phaseLeft, *phaseRight, *fifoOutput;
//float sum_ensig;
//float beta;
int on;
int frameCounter=0;
float snr_counter=0;
float sum_SNR=0;
float snr_avg;
float musical=2;
int numchan = 9;
//int frameCounter=0;
float chanbuf[9][512];
float *chan_rms = (float*)calloc(numchan,sizeof(float));
float *avg_rms = (float*)calloc(numchan,sizeof(float));
float *gain = (float*)calloc(numchan,sizeof(float));
float *gain_L = (float*)calloc(numchan,sizeof(float));
float *outputR = (float*)calloc(nFFT,sizeof(float));
float *outputL = (float*)calloc(nFFT,sizeof(float));


@implementation ViewController {
        SuperpoweredIOSAudioIO *audioIO;
    SuperpoweredFrequencyDomain *frequencyDomain;
   // float *output;
    float *magnitudeLeft, *magnitudeRight, *phaseLeft, *phaseRight, *fifoOutput;
        int fifoOutputFirstSample, fifoOutputLastSample, stepSize, fifoCapacity;
    
    
}

@synthesize betaLabel;
@synthesize betaSlider;
@synthesize snrLabel;
@synthesize EnhancedSwitch;
@synthesize OnDisplay;
@synthesize muuLabel;
@synthesize vLabel;
@synthesize musicalLabel;
@synthesize muuSlider;
@synthesize vSlider;
@synthesize musicalStepper;
@synthesize musical;
@synthesize muu;
@synthesize v;
@synthesize stepper;

-(void) updateBeta{
    x=stepper.value;
    betaLabel.text=[ NSString stringWithFormat:@"%f",x];
    }
/*-(void) updatemuu{
    x1=muuSlider.value;
    muu.text=[ NSString stringWithFormat:@"%f",x1];
}
-(void) updatev{
    x2=vSlider.value;
    v.text=[ NSString stringWithFormat:@"%f",x2];
}
-(void) updatemusical{
    x3=musicalStepper.value;
    musical.text=[ NSString stringWithFormat:@"%f",x3];
}*/


- (IBAction)betaValue:(id)sender {
   
    }
/*- (IBAction)musicalValue:(id)sender {
    [self updatemusical];
}

- (IBAction)muuValue:(id)sender {
    [self updatemuu];
}

- (IBAction)vValue:(id)sender {
    [self updatev];
}*/


- (IBAction)switchPressed:(id)sender {
    if(EnhancedSwitch.on){
        on=1;
        OnDisplay.text=[ NSString stringWithFormat:@"ON"];
    }
        else{
            on=0;
            OnDisplay.text=[ NSString stringWithFormat:@"OFF"];
            
    }
}
float rms(float *input,int size)
{
    float total = 0;
    for(int i=0;i<size;i++)
        if(input[i]<20&&input[i]>0) {
            total += input[i] * input[i];
        }
    
    return total/size;
}

float compute_gain(float rms_temp, float TK_temp){
    float current_dB,adjust_dB;
    float gain_temp ;
    current_dB = 20 * log10(rms_temp / 0.00002);
    adjust_dB = TK_temp - current_dB;
    gain_temp = pow(10, 0.05 * adjust_dB);
    gain_temp = 1.1;
    return gain_temp;
}

#define FFT_LOG_SIZE 10 // 2^11 = 2048
// This callback is called periodically by the audio system.
static bool audioProcessing(void *clientdata, float **buffers, unsigned int inputChannels, unsigned int outputChannels, unsigned int numberOfSamples, unsigned int samplerate, uint64_t hostTime)
{
    __unsafe_unretained ViewController *self = (__bridge ViewController *)clientdata;
        // Input goes to the frequency domain.
    float interleaved[numberOfSamples * 2 + 16];
    SuperpoweredInterleave(buffers[0], buffers[1], interleaved, numberOfSamples);
    self->frequencyDomain->addInput(interleaved, numberOfSamples);
    
       // In the frequency domain we are working with 1024 magnitudes and phases for every channel (left, right), if the fft size is 2048.
    while (self->frequencyDomain->timeDomainToFrequencyDomain(self->magnitudeLeft, self->magnitudeRight, self->phaseLeft, self->phaseRight))
    {
    beta=self->x;
      //  muu=self->x1;
      //  v=self->x2;
      //  musical=self->x3;
        //printf("%f\n",v);
        float total_noisepower=0;
        float total_speechpower=0;
        float sum_ensig = 0;
        float sum_nsig = 0;
        if(on==1)
        {
            if(beta==0)
            {
                beta=0.5;
            }
          
            frameCounter++;
            if(frameCounter<=7)
            {
                for (int i = 0; i < nFFT; i++)
                {
                    noise_mean[i] += self->magnitudeLeft[i];
                   // printf("%f\n",noise_mean[i]);
                }
            }
            if(frameCounter==7)
            {
                for (int i = 0; i < nFFT; i++)
                {
                    noise_mu2[i]=(float)pow(noise_mean[i]/6,2);
                   // printf("%f\n",noise_mu2[i]);
                }
            }
            if(frameCounter>=7)
            {
                for (int i = 0; i < nFFT; i++)
                    sig2[i] = (float)pow(self->magnitudeLeft[i], 2);
                
                for (int i = 0; i < nFFT; i++)
                {
                    if( noise_mu2[i]==0)
                        noise_mu2[i] = noise_mu2[i] + epsilon;
                    gammak[i] = sig2[i] / noise_mu2[i] < 40 ? sig2[i] / noise_mu2[i] : 40;
                   // printf("%f\n",gammak[i]);
                    
                }
                if (frameCounter == 7)
                {
                    for (int i = 0; i < nFFT; i++) {
                        max = gammak[i] - 1 > 0 ? gammak[i] - 1 : 0;
                        ksi[i] = aa + (1 - aa) * max;
                       // printf("%f\n",ksi[i]);

                    }
                    
                }
                else
                {
                    for (int i = 0; i < nFFT; i++)
                    {
                        if( noise_mu2[i]==0)
                            noise_mu2[i] = noise_mu2[i] + epsilon;
                        max = gammak[i] - 1 > 0 ? gammak[i] - 1 : 0;
                        ksi[i] = aa * (Xk_prev[i] / noise_mu2[i]) + (1 - aa) * max;
                        ksi[i] = ksi_min > ksi[i] ? ksi_min : ksi[i];
                       // printf("%f\n",Xk_prev[i]);
                    }
                }
               
                sum_log_sigma_k = 0;
                for (int i = 0; i < nFFT; i++)
                {
                    log_sigma_k[i] = gammak[i] * ksi[i] / (1 + ksi[i]) - log(1 + ksi[i]);
                    if(log_sigma_k[i]<1000)
                    {
                        sum_log_sigma_k += log_sigma_k[i];
                    }
                    
                
                }
               // printf("%f\n",sum_log_sigma_k);
                vad_decision = sum_log_sigma_k / (nFFT);
                
                //__android_log_print(ANDROID_LOG_VERBOSE, "MyApp", "The value of SNR_db is %f",vad_decision );
                if (vad_decision < eta)// % noise on
                {
                    
                    for (int i = 0; i < nFFT; i++)
                    {
                        noise_pow[i] = noise_pow[i] + noise_mu2[i];
                    }
                    
                    count = count + 1;
                    for (int i = 0; i < nFFT; i++)
                        noise_mu2[i] = noise_pow[i] / count;
                }
                
                
                for (int i = 0; i < nFFT; i++)
                    vk[i] = ksi[i] * gammak[i] / (1 + ksi[i]);
                
                for (int i = 0; i < nFFT; i++)
                {
                    hw[i] = (ksi[i] + sqrt(pow(ksi[i], 2) + ((2*beta)*(beta + ksi[i])) * ksi[i] / (gammak[i]+epsilon))) / (2*beta * (beta + ksi[i]));

                  // printf("%f\n",hw[i]);
                    
                   // hw[i] = (ksi[i] + (float)sqrt(pow(ksi[i], 2) + (1 + ksi[i] / beta) * ksi[i] / (gammak[i]+epsilon))) / (2 * (beta + ksi[i]));
                }
            
               for (int i = 0; i < nFFT; i++)
                {
                    if(ensig[i]<100)
                    {
                    sum_ensig += pow(ensig[i], 2);
                    }
                    if(self->magnitudeLeft[i]<100)
                    {
                    sum_nsig += pow(self->magnitudeLeft[i], 2);
                    }
                }
              //  printf("%f\n",sum_ensig);
              //  printf("%f\n",sum_nsig);

           float PR = sum_ensig/sum_nsig;
                
              if (PR >= 0.4)
                    PRT = 1;
                else
                    PRT = PR;
                
                if (PRT == 1)
                    N = 1;
                else
                    N = 2 * round((1 - PRT / 0.4) * 5) + 1;
                
                //printf("%f\n",N);
                for (int i = 0; i < (int) N; i++)
                    H[i] = 1 / N;
                //	H(1:N) = 1 / N;
                
                for (int i = 0; i < N + nFFT - 1; i++)
                {
                    int kmin, kmax, k;
                    
                    HPF[i] = 0;
                    
                    kmin = (i >= nFFT - 1) ? i - (nFFT - 1) : 0;
                    kmax = (i < N - 1) ? i : N - 1;
                    
                    for (k = kmin; k <= kmax; k++)
                        HPF[i] += H[k] * sqrt(hw[i - k]* hw[i - k]);
                    //printf("%f\n",HPF[i]);
                    //hw[i]=HPF[i];
                    
                   
                    
                    
                
                if(beta<=0.8)
                {
                    HPF[i]=1;
                }
                }
                if (SPU == 1)
                    
                    for (int i = 0; i < nFFT; i++)
                    {
                        evk[i] = exp(vk[i]);
                        Lambda[i] = qkr * evk[i] / (1 + ksi[i]);
                        pSAP[i] = Lambda[i] / (1 + Lambda[i]);
                        self->magnitudeLeft[i] = self->magnitudeLeft[i] * HPF[i] * pSAP[i];
                    }
                else
                    for (int i = 0; i < nFFT; i++)
                    {
                       self->magnitudeLeft[i] = self->magnitudeLeft[i] * HPF[i];
                       ensig[i]=self->magnitudeLeft[i];
                       self->magnitudeRight[i] = self->magnitudeLeft[i];
                        if (self->magnitudeLeft[i]<10)
                        {
                            // total_speechpower+=magnitudeLeft[i]*magnitudeLeft[i];
                            total_speechpower+=(float)pow(self->magnitudeLeft[i],2);
                        }
                        if(noise_pow[i]<100)
                        {
                            total_noisepower+=noise_pow[i];
                        }
                    }
                
                for (int i = 0; i < nFFT; i++)
                {
                    Xk_prev[i] =(float) pow(self->magnitudeLeft[i],2);
                    //X->real[i]=X->real[i]*HPF[i];
                    //X->imaginary[i]=X->imaginary[i]*HPF[i];
               // printf("%f\n",hw[i]);
                }

            }
            SNR_db=20*log((2* total_speechpower)/total_noisepower);

            /*CALCULATE AVERAGE OF SNR*/
            
            sum_SNR+=SNR_db;
            snr_counter++;
            if (snr_counter==100)
            {
                //pthread_mutex_lock(&mutex);
                snr_avg=sum_SNR/snr_counter;
                //  __android_log_print(ANDROID_LOG_VERBOSE, "MyApp", "The value of SNR_db is %f",snr_counter );
                sum_SNR=0;
                snr_counter=0;
            }
            
        ///////****************** SE ends here; output is self->magnitudeLeft and self->magnitudeRight *************//////
            
            
            

            for (int j = 0; j < numchan; j++) {
                for(int i=0;i<nFFT;i++){
                    chanbuf[j][i] = self->magnitudeLeft[i] * FB[j][i];
                }
                /////
                avg_rms[j] = rms(chanbuf[j], nFFT)*0.2 + avg_rms[j]*0.8;
                gain[j] = compute_gain(avg_rms[j],TK[j]);
                gain_L[j] = compute_gain(avg_rms[j],TK_L[j]);
                                /*
                 current_dB = 20*log10(avg_rms[j]/0.00002);
                 adjust_dB = TK[j] - current_dB;
                 gain[j] = pow(10, 0.05*adjust_dB);
                 adjust_dB_L = TK_L[j] - current_dB;
                 gain_L[j] = pow(10, 0.05*adjust_dB_L);
                 if (C_flag == 1)
                 {
                 gain[j] = 1.15;
                 gain_L[j] = 1.15;
                 }
                 */
                //////
                for(int i=0;i<nFFT;i++){
                    if(j==0) {
                        outputR[i] = chanbuf[j][i] * gain[j];
                        outputL[i] = chanbuf[j][i] * gain_L[j];
                    }
                    else {
                        outputR[i] += chanbuf[j][i]*gain[j];
                        outputL[i] += chanbuf[j][i]*gain_L[j];
                    }
                }
                
            }
            for(int i=0;i<nFFT;i++){
                self->magnitudeRight[i] = outputR[i] * gainR[i];
                self->magnitudeLeft[i] = outputL[i] * gainL[i];
//                self->magnitudeRight[i] = self->magnitudeRight[i];
//                self->magnitudeLeft[i] = self->magnitudeLeft[i] ;
                
            }
       
            
            
            
            
            
            
            
            
            
            
    
        }
    
   // return true;
               // pthread_mutex_unlock(&mutex);
                // return true;
        
            // This is just a quick example: we remove the magnitude of the first 20 bins, meaning total bass cut between 0-430 Hz.
      //memset(self->magnitudeLeft, 0, 80);
        //memset(self->magnitudeRight, 0, 80);

        // We are done working with frequency domain data. Let's go back to the time domain.

        // Check if we have enough room in the fifo buffer for the output. If not, move the existing audio data back to the buffer's beginning.
       if (self->fifoOutputLastSample + self->stepSize >= self->fifoCapacity)  // This will be true for every 100th iteration only, so we save precious memory bandwidth.
        {
            int samplesInFifo = self->fifoOutputLastSample - self->fifoOutputFirstSample;
            if (samplesInFifo > 0) memmove(self->fifoOutput, self->fifoOutput + self->fifoOutputFirstSample * 2, samplesInFifo * sizeof(float) * 2);
            self->fifoOutputFirstSample = 0;
            self->fifoOutputLastSample = samplesInFifo;
        };

        // Transforming back to the time domain.
        self->frequencyDomain->frequencyDomainToTimeDomain(self->magnitudeLeft, self->magnitudeRight, self->phaseLeft, self->phaseRight, self->fifoOutput + self->fifoOutputLastSample * 2);
        self->frequencyDomain->advance();
        self->fifoOutputLastSample += self->stepSize;
    
    };

    // If we have enough samples in the fifo output buffer, pass them to the audio output.
    if (self->fifoOutputLastSample - self->fifoOutputFirstSample >= numberOfSamples) {
        SuperpoweredDeInterleave(self->fifoOutput + self->fifoOutputFirstSample * 2, buffers[0], buffers[1], numberOfSamples);
        // buffers[0] and buffer[1] now have time domain audio output (left and right channels)
        self->fifoOutputFirstSample += numberOfSamples;
        return true;
    } else return false;
}

- (void)viewDidLoad {
    [super viewDidLoad];
    OnDisplay.text=@"OFF";
betaLabel.text = @"0.5";
    snrLabel.text=@"0";
    muu.text = @"0.5";
    v.text = @"0.01";
    musical.text = @"2";
    frequencyDomain = new SuperpoweredFrequencyDomain(FFT_LOG_SIZE); // This will do the main "magic".
    stepSize = frequencyDomain->fftSize / 4; // The default overlap ratio is 4:1, so we will receive this amount of samples from the frequency domain in one step.

    // Frequency domain data goes into these buffers:
    magnitudeLeft = (float *)malloc(frequencyDomain->fftSize * sizeof(float));
    magnitudeRight = (float *)malloc(frequencyDomain->fftSize * sizeof(float));
    phaseLeft = (float *)malloc(frequencyDomain->fftSize * sizeof(float));
    phaseRight = (float *)malloc(frequencyDomain->fftSize * sizeof(float));

    // Time domain result goes into a FIFO (first-in, first-out) buffer
    fifoOutputFirstSample = fifoOutputLastSample = 0;
    fifoCapacity = stepSize * 100; // Let's make the fifo's size 100 times more than the step size, so we save memory bandwidth.
    fifoOutput = (float *)malloc(fifoCapacity * sizeof(float) * 2 + 128);

    // Audio input/output handling.
    audioIO = [[SuperpoweredIOSAudioIO alloc] initWithDelegate:(id<SuperpoweredIOSAudioIODelegate>)self preferredBufferSize:12 preferredMinimumSamplerate:44100 audioSessionCategory:AVAudioSessionCategoryPlayAndRecord channels:2 audioProcessingCallback:audioProcessing clientdata:(__bridge void *)self];
    [audioIO start];
}

- (void)dealloc {
    audioIO = nil;
    delete frequencyDomain;
    free(magnitudeLeft);
    free(magnitudeRight);
    free(phaseLeft);
    free(phaseRight);
    free(fifoOutput);
}

- (void)interruptionStarted {}
- (void)interruptionEnded {}
- (void)recordPermissionRefused {}
- (void)mapChannels:(multiOutputChannelMap *)outputMap inputMap:(multiInputChannelMap *)inputMap externalAudioDeviceName:(NSString *)externalAudioDeviceName outputsAndInputs:(NSString *)outputsAndInputs {}
- (IBAction)snrDisplay:(id)sender{
    
    snrLabel.text=[NSString stringWithFormat:@"%f",SNR_db];

}
- (IBAction)stepbeta:(id)sender {
    [self updateBeta];
}
@end

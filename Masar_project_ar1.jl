#Creating AR1-process randomly
#using Pkg
#Pkg.add("Distributions")
using Random, Distributions
function AR_1(alpha,sigma2,iter)
    y=[0.0 for i=1:iter]
    d=Normal(0,sqrt(sigma2))
    epsilon=rand(d,iter)
    for i in 1:iter-1
        y[i+1]=alpha*y[i]+epsilon[i]
    end
    return y
end
n=100
sigma2=1
alpha_start=0.5
AR_y=AR_1(alpha_start,sigma2,n)
#Pkg.add("Plots")
using Plots
xaxis=1:n
plot(xaxis,AR_y, title="AR(1)-process", label="AR(1)", xlabel="time")

#MLE for AR(1)
function MLE(x)
    alpha_MLE_s1=0
    alpha_MLE_s2=0
    mean_ar=mean(x)
    for i in 1:n-1
        alpha_MLE_s1=alpha_MLE_s1+(x[i]-mean_ar)*(x[i+1]-mean_ar)
        alpha_MLE_s2=alpha_MLE_s2+(x[i]-mean_ar)^2
    end
    return alpha_MLE_s1/alpha_MLE_s2
end
alpha_MLE=MLE(AR_y)
#Log-likelihood function for the time series
function log_likelihood(sigma2,alpha)
    x=AR_y
    n=length(x)
    temp_cond=-(n-1)/(2)*log(2*pi)-(n-1)/(2)*log(sigma2)
    temp_marg=-(n/2)*log(2*pi)-(1/2)*log(sigma2/(1-alpha^2))-(1-alpha^2)/(2*sigma2)*x[1]
    temp_faktor_cond=-1/(2*sigma2)
    temp_sum_cond=0
    for i in 2:n
        temp_sum_cond=temp_sum_cond+(x[i]-alpha*x[i-1])^2
    end
    return temp_cond+temp_marg+temp_faktor_cond*temp_sum_cond
end
alpha_values=0:0.9/(n/2):0.9
y_alpha=log_likelihood.(sigma2,alpha_values)
plot(alpha_values,y_alpha, title="Log-likelihood", label="Log-likelihood", xlabel="Time domain")
alpha_values[findmax(y_alpha)[2]]

#Spectrum
using FFTW
function AR_1_spectrum(omega, alpha, sigma2)
    #y=sigma2*(1-alpha^2)/(pi*(1-2*alpha*cos(omega)+alpha^2))
    y=sigma2/(pi*(1-2*alpha*cos(omega)+alpha^2))
    return y
end
omega=rfftfreq(n)*2*pi
y_freq=AR_1_spectrum.(omega,alpha_start,sigma2)
plot(omega,y_freq, title="Spectrum of AR(1)", label="Spectral density function", xlabel="Frequency domain")

#FFTW
#using Pkg
#Pkg.add("FFTW")
using FFTW
function fast_fourier_transform(x)
    y_freq_fft=rfft(x)
    y_freq_fft_real=abs.(y_freq_fft).^(2)./(2*pi*n)#abs.(y_freq_fft).^(2)./(2*pi*n)
    return y_freq_fft_real
end
omega_fft=rfftfreq(n)*2*pi
y_freq_fft_real=fast_fourier_transform(AR_y)
plot(omega_fft,y_freq_fft_real, title="Periodogram", label="Periodogram", xlabel="Frequency domain")
plot!(omega_fft,y_freq, label="Spectral density function", xlabel="Frequency domain")

#Whittle maximum likelihood, with FFT
alpha_values_fft=alpha_values
function log_whittle_like_fft(alpha)
    n=length(AR_y)
    omega=rfftfreq(n)*2*pi
    s_sum=0
    for i in 1:length(omega)
        f_omega=AR_1_spectrum(omega[i], alpha, sigma2)
        #f_omega=sigma2*(1-alpha^2)/(pi*(1-2*alpha*cos(omega[i])+alpha^2))
        I_omega=y_freq_fft_real[i]
        s_sum=s_sum+log(f_omega)+I_omega/f_omega
    end
    w=(-1/2)*(2*n*log(2*pi)+2*s_sum)
    return w
end
whittle_like_fft=log_whittle_like_fft.(alpha_values_fft)
plot(alpha_values_fft,whittle_like_fft,title="Log-likelihood", label="Whittle log-likelihood", xlabel="Alpha")
whittle_MLE_fft=alpha_values_fft[argmax(whittle_like_fft)]
#Normalize likelihood
plot(alpha_values_fft,exp.(whittle_like_fft.-findmax(whittle_like_fft)[1]),title="Log-likelihood", label="Whittle log-likelihood", xlabel="Alpha")
plot!(alpha_values,exp.(y_alpha.-findmax(y_alpha)[1]), title="Log-likelihood", label="Log-likelihood", xlabel="Alpha")
#Save plot for reporting
savefig("/Users/masaral-mosawi/Documents/VS Code/Julia plots/Whittle_test.pdf")
#Credible interval periodogram, simulated
function cred_int(n)
    c=zeros(length(omega_fft))
    for i in 1:n
        c=hcat(c,fast_fourier_transform(AR_1(alpha_start,sigma2,n)))
    end
    return c
end
cred_interval=cred_int.(n+1)
q_upper=zeros(length(omega_fft))
q_lower=zeros(length(omega_fft))
for i in 1:length(omega_fft)
    q_lower[i]=quantile(cred_interval[i,:], [0.025,0.975])[1]
    q_upper[i]=quantile(cred_interval[i,:], [0.025,0.975])[2]
end
plot(omega_fft,y_freq_fft_real, title="Periodogram using different methods", label="FFT", xlabel="Frequency domain")
plot!(omega,y_freq, label="Spectral density function", xlabel="Frequency domain")
plot!(omega,q_lower, label="Lower sim credible interval", xlabel="Frequency domain")
plot!(omega,q_upper, label="Upper sim credible interval", xlabel="Frequency domain")
#Credible interval, P(Chisq2,1-alpha/2<2I(w)/f(w)<Chisq2,alpha/2)=1-alpha
chi2_lower=(quantile(Chisq(2), 0.025))
chi2_upper=(quantile(Chisq(2), 0.975))
plot(omega_fft,y_freq_fft_real, title="Periodogram using different methods", label="FFT", xlabel="Frequency domain")
plot!(omega,y_freq, label="Spectral density function", xlabel="Frequency domain")
plot!(omega_fft,y_freq.*chi2_lower/2, label="Lower theor credible interval", xlabel="Frequency domain")
plot!(omega_fft,y_freq.*chi2_upper/2, label="Upper theor credible interval", xlabel="Frequency domain")
plot!(omega,q_upper, label="Upper sim credible interval", xlabel="Frequency domain")
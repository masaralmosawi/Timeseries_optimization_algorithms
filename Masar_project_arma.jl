#Creating ARIMA observations, https://github.com/JuliaDynamics/ARFIMA.jl
# Masar 202102: ] add https://github.com/JuliaDynamics/ARFIMA.jl, in Julia console window
using ARFIMA
n=1000
sigma2=1
#p=promote(0.273, -0.81)#AR(p), Check: so the arma(1,1) is correct against known plots
#q=promote(-0.9)#MA(q)
p=promote(0.5)#AR(p)
q=promote(0.5)#MA(q)
ARMA_y=arfima(n,sqrt(sigma2),nothing,SVector(p), SVector(q))
#Pkg.add("Plots")
using Plots
xaxis=1:n
plot(xaxis,ARMA_y, title="ARMA(p,q)-process", label="ARMA", xlabel="Time")

#Spectrum
using FFTW
function ARMA_spectrum(omega, sigma2)
    #B_phi_vector=[exp(-2*pi*omega*(1im)*i) for i=1:length(p)]
    #B_theta_vector=[exp(-2*pi*omega*(1im)*i) for i=1:length(q)]
    B_phi_vector=[exp(-omega*(1im)*i) for i=1:length(p)]
    B_theta_vector=[exp(-omega*(1im)*i) for i=1:length(q)]
    phi_B=p.*B_phi_vector
    theta_B=q.*B_theta_vector
    y=sigma2*(abs(1-sum(theta_B))^2)/(abs(1-sum(phi_B))^2)
    return y
end
omega=rfftfreq(n)*2*pi
y_freq=ARMA_spectrum.(omega,sigma2)
plot(omega,y_freq, title="Spectrum of ARMA(p,q)", label="Spectral density function", xlabel="Frequency domain")

#FFTW
#using Pkg
#Pkg.add("FFTW")
using FFTW
function fast_fourier_transform(x)
    y_freq_fft=rfft(x)
    y_freq_fft_real=abs.(y_freq_fft).^(2)./(2*pi*n)
    return y_freq_fft_real
end
y_freq_fft_real=fast_fourier_transform(ARMA_y)
plot(omega,y_freq_fft_real, title="Periodogram", label="FFT", xlabel="Frequency domain")
plot!(omega,y_freq, label="Spectral density function", xlabel="Frequency domain")

#Contour Plots
using FFTW
function ARMA_spectrum_function(omega, sigma2, p, q)
    #B_phi_vector=[exp(-2*pi*omega*(1im)*i) for i=1:length(p)]
    #B_theta_vector=[exp(-2*pi*omega*(1im)*i) for i=1:length(q)]
    B_phi_vector=[exp(-omega*(1im)*i) for i=1:length(p)]
    B_theta_vector=[exp(-omega*(1im)*i) for i=1:length(q)]
    phi_B=p.*B_phi_vector
    theta_B=q.*B_theta_vector
    y=sigma2*(abs(1-sum(theta_B))^2)/(abs(1-sum(phi_B))^2)
    return y
end
#y_freq_function=ARMA_spectrum_function.(omega,sigma2, 0.2,0)
#plot(omega,y_freq_function, title="Spectrum of ARMA(p,q)", label="Spectral density function", xlabel="Frequency domain")
#Whittle maximum likelihood, function
function log_whittle_like_fft_function(phi, theta)
    n=length(ARMA_y)
    omega=rfftfreq(n)*2*pi
    s_sum=0
    for i in 1:length(omega)
        f_omega=ARMA_spectrum_function(omega[i], sigma2, phi, theta)
        I_omega=y_freq_fft_real[i]
        s_sum=s_sum+log(f_omega)+I_omega/f_omega
    end
    w=(-1/2)*(2*n*log(2*pi)+2*s_sum)
    return w
end
a=-0.7
b=0.7
#n=length(y_freq_function)
phi_values_func=a:(b-a)/(n-1):b
theta_values_func=a:(b-a)/(n-1):b
z=zeros(n,n)
for i in 1:n
    for j in 1:n
        z[i,j]=log_whittle_like_fft_function(phi_values_func[i], theta_values_func[j])
        #z[i,j]=phi_values_func[i]^2-theta_values_func[j]^2
    end
end
#contour(phi_values_func,theta_values_func,z.-findmax(z)[1],clim=(-10,0), title="Contour of ARMA(1,1)", xlabel="Phi", ylabel="Theta")
contour(theta_values_func,phi_values_func,z.-findmax(z)[1],clim=(-10,0), title="Contour of ARMA(1,1)", xlabel="Theta(q)", ylabel="Phi(p)")
#Save plot for reporting
savefig("/Users/masaral-mosawi/Documents/VS Code/Julia plots/contour_root_canc.pdf")

#Likelihood with R

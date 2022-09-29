import numpy as np
import timeit

time_start = timeit.default_timer() #starter tiden


m = 0.169 #definerer massen til kula i kg
R = 0.50 #definerer radien R til kvartsirkelen i m
r = 0.057/2 #definerer radien r til kula i m

start_omega = 0 #definerer vinkelfarten ved t=0
start_alpha = 0 #definerer vinkelakselerasjonen ved t=0
start_theta = 4*np.pi/180   #definerer startvinkelen ved t=0
delta_t = 0.00001 #kan ikke være 0 for å kunne bruke Eulers metode


def angle(m, R, r, omega_prior, alpha, theta_prior, delta_t):
    """
    Funksjon som tar inn masse, en gitt radius R til kvartsirkelen, radiusen r til et objekt, 
    en vinkelfart, vinkelakselerasjon, startvinkel og tidsintervall delta_t.
    
    Deretter returnerer funksjonen en streng "func_values" med informasjon om hvilken vinkel
    objektet forlater banen i. Dette kommer i både grader og radianer.
    Vi får også vite antall iterasjoner, samt tiden ballen brukte på å forlate banen.
    Vi får også energiendringen fra start til slutt i joule og prosent.
    
    """
    g = 9.81 #definerer tyngdeakselerasjonen i m/s^2
  
    energy_prior = m*g*np.cos(theta_prior)*(R+r) #startenergi justert for startvinkel

    N = 1 #definerer en vilkårlig verdi for N større enn null slik at løkken vil kjøre
    i = 0 #setter i for å telle antall iterasjoner  
    while N > 0: #while løkke som itererer frem til N blir mindre enn 0. 
        alpha = g * np.sin(theta_prior) / (r + R) #regner ut vinkelakselerasjonen
        omega_post = omega_prior + alpha * delta_t #bruker Eulers metode for å finne  vinkelfarten
        theta_post = theta_prior + omega_prior * delta_t #bruker Eulers metode for å finne vinkelen
        
        v = (r + R) * omega_prior #farten langs overflaten
        N = m * g * np.cos(theta_prior) - (m * (v ** 2)) / (r + R) #regner ut normalkraften
        
        omega_prior = omega_post #oppdaterer omega før neste iterasjon
        theta_prior = theta_post #oppdater theta før neste iterasjon
        
        i += 1 #legger til en i iterasjonstellingen
        
    pot_e =  m*g*(np.cos(theta_post)*(R+r)) #regner ut potensiell energi ved slutt
    kin_e = 0.5 * m * v ** 2 #regner ut kinetisk energi ved slutt
        
    energy_post = pot_e + kin_e #regner ut den totalte energien ved slutten
        
    func_values = f"""
Objektet forlot banen etter {theta_post*180/np.pi} grader,
som er {theta_post} radianer. 

Løkken kjørte {i} ganger med delta_t: {delta_t}
Objektet forlot banen etter {round(delta_t*i,6)} sekunder.

Energiendringeen ble {(energy_post - energy_prior)} joule, dette tilsvarer 
en endring i energi på {(energy_post - energy_prior) / energy_prior * 100}%
"""#lager en stor streng så svaret blir printet fint
    return func_values #returnerer strengen med verdiene

    
print(angle(m,R,r,start_omega,start_alpha, start_theta, delta_t)) #kjører og printer funksjonen

time_stop = timeit.default_timer() #stopper tiden

print(f"Progammet brukte {time_stop-time_start} sekunder") #printer tiden programmet brukte for å kjøre
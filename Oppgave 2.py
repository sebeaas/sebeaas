import numpy as np
import matplotlib.pyplot as plt
import timeit

time_start = timeit.default_timer() #starter tiden

g = 9.81 #definerer tyngdeakselerasjonen i m/s^2
R = 0.50 #definerer radien R til kvartsirkelen i m

omega_prior = 0 #definerer vinkelhastigheten ved t=0 
alpha = 0 #definerer vinkelakselerasjonen ved t=0
theta_prior = 6.6*np.pi/180  #definerer startvinkelen ved t=0
delta_t = 0.00001 #kan ikke være 0 for å kunne bruke Eulers metode

objects = {'kule': [2/5, 0.169, 0.057/2], 'kompakt sylinder': [1/2, 1.099, 0.044/2]} 
#dictionary for objektene vi bruker i denne oppgaven


def angle(objects, object, R, omega_prior, alpha, theta_prior):
    
    """
    Funksjon som tar inn masse, en gitt radius R til kvartsirkelen, radiusen r til et objekt, 
    en vinkelfart, vinkelakselerasjon, startvinkel og tidsintervall delta_t.
    
    Deretter returnerer funksjonen en streng "func_values" med informasjon om hvilken vinkel
    objektet forlater banen i. Dette kommer i både grader og radianer.
    Vi får også vite antall iterasjoner, samt tiden ballen brukte på å forlate banen.
    Vi får også energiendringen fra start til slutt i joule og i prosent.
    
    Som tilleggsinformasjon for å verifisere at løkken kjører riktig får vi et plot som 
    inneholder objektets bane med x,y koordinater. 
    """
        
    values = objects.get(object) #lagrer informasjon om et gitt objekt i en variabel
    c = values[0] #henter ut c fra dict
    m = values[1] #henter ut m fra dict
    r = values[2] #henter ut r fra dict
        
    
    x_0 = 0 #definerer startverdi for x 
    y_0 = (R+r)*np.cos(theta_prior) #definerer startverdi for y, tar hensyn til startvinkel
    v = 0 #definerer startfarten i t=0
    
    pos_arr_y = np.zeros(2*int(1/delta_t)) 
    pos_arr_x = np.zeros(2*int(1/delta_t)) 
    #definerer to zero-arrays, et for x-posisjon og et for y-posisjon, legger også til en god
    #sikkerhetsmargin slik at løkken ikke vil bryte av at indeksen ikke strekker til. 
    #Forhåndsdefinerer arrayene for at løkken skal kjøre raskt og bruke mindre minne.

    
    energy_prior = m*g*np.cos(theta_prior)*(R+r) #startenergi med høyde justert for startvinkel
        
    N = 1 #definerer en vilkårlig verdi for N større enn null slik at løkken vil kjøre
    i = 0 #setter i for å telle antall iterasjoner  
    
    while N > 0: #løkke som kjører til normalkraften er mindre enn 0
        alpha = ((g * np.sin(theta_prior)) / (1 + c)) / (r + R) #vinkelakselerasjon
        omega_post = omega_prior + alpha * delta_t #bruker Eulers metode for å finne vinkelfart
        theta_post = theta_prior + omega_prior * delta_t #bruker Eulers metode for å finne vinkelen
        
        x = x_0 + (((R+r) * omega_prior * np.cos(theta_post) ) + v * np.cos(theta_prior)) * delta_t / 2 #regner ut x-posisjonen ved hjelp av S=1/2(v-v_0)*t, bruker v fra forrige iterasjon som v_0, og regner ut nye v direkte
        pos_arr_x[i] = x #legger til x-posisjonen i x-arrayet
        x_0 = x #oppdaterer x posisjonen
        
        y = y_0 - ((R+r) * omega_prior * np.sin(theta_post) + v * np.sin(theta_prior)) * delta_t / 2 #regner ut y-posisjonen ved hjelp av S=1/2(v-v_0)*t, bruker v fra forrige iterasjon som v_0, og regner ut nye v direkte
        pos_arr_y[i] = y #legger til x-posisjonen i x-arrayet
        y_0 = y #oppdaterer y-posisjonen
        

        v = (r + R) * omega_prior #farten langs overflaten
        N = m * g * np.cos(theta_prior) - (m * (v ** 2)) / (r + R) #regner ut normalkraften
        
        i += 1 #legger til en i antall iterasjoner
        
        omega_prior = omega_post #oppdaterer omega 
        theta_prior = theta_post #oppdater theta 
        
    pot_e =  m*g*(np.cos(theta_post)*(R+r)) #regner ut potensiell energi ved slutten
    kin_e = 0.5 * m * v ** 2 #regner ut kinetisk energi ved slutten
    
    I = c * m * r ** 2 #regner ut dreiemomentet

    rot_e = 0.5 * I * (v/r) ** 2 #regner ut rotasjonsenergien
        
    energy_post = pot_e + kin_e + rot_e #regner ut total energi
    
    plot_x = np.delete(pos_arr_x, [*range(i,2*int(1/delta_t))]) #lager nytt array uten nullene fra arrayet som ikke tas i bruk
    plot_y = np.delete(pos_arr_y, [*range(i,2*int(1/delta_t))]) #lager nytt array uten nullene fra arrayet som ikke tas i bruk
    
    plt.title(f"Banen til {object}") #plotter tittel
    plt.xlabel("x-posisjon") #plotter navn på x-aksen
    plt.ylabel("y-posisjon") #plotter navn på y-aksen
    plt.plot(plot_x, plot_y, "r-") #plotter x og y arrayene
        
    plt.grid() #legger til grid
    plt.show() #viser plottet
    
    func_values = f"""
Ballen forlot banen etter {theta_post*180/np.pi} grader, som er {theta_post} radianer.
Løkken kjørte {i} ganger med delta_t: {delta_t}
Ballen forlot banen etter {round(delta_t*i,6)} sekunder.
Energiendringeen ble {(energy_post - energy_prior)} joule, dette tilsvarer en 
endring i energi på {(energy_post - energy_prior) / energy_prior * 100}%
Som tilleggsinformasjon for å verifisere at løkken kjører riktig får vi et plot som 
inneholder objektets bevegelse med x,y koordinater.
"""#lager en stor streng så svaret blir printet fint
    return func_values
    #returnerer vinkelen i grader, samt antall kjøringer, normalkraften
    #der objektet mister kontakt med underlaget, og den siste vinkelen registrert

print(angle(objects, "kompakt sylinder", R, omega_prior, alpha, theta_prior)) #kjører og printer funksjonen

time_stop = timeit.default_timer() #stopper tiden

print(f"Progammet brukte {time_stop-time_start} sekunder") #printer tiden programmet brukte for å kjøre
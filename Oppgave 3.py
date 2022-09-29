import numpy as np
import matplotlib.pyplot as plt
import timeit

time_start = timeit.default_timer() #starter tiden

objects = {'hul sylinder': [1.0, 0.254, [0.042/2, 0.036/2]], 'kompakt sylinder': [1/2, 1.099, [0.044/2]]} 
#dictionary for objektene vi bruker i denne oppgaven

R = 0.5 #radien til kvartsirkelen
omega_prior = 0 #definerer vinkelfarten ved t=0
alpha = 0 #definerer vinkelakselerasjonen ved t=0
theta_prior = 4.1 * np.pi / 180 #definerer startvinkelen ved t=0
delta_t = 0.00001 #kan ikke være 0 for å kunne bruke Eulers metode

def angle(objects, object, R, omega_prior, alpha, theta_prior, delta_t):
    
    """
    Funksjon som tar inn masse, en gitt radius R til kvartsirkelen, radiusen r til et objekt, 
    en vinkelfart, vinkelakselerasjon, startvinkel og tidsintervall delta_t.
    
    Deretter returnerer funksjonen en streng "func_values" med informasjon om hvilken vinkel
    objektet forlater banen i. Dette kommer i både grader og radianer.
    Vi får også vite antall iterasjoner, samt tiden ballen brukte på å forlate banen.
    Vi får også energiendringen fra start til slutt i joule og i prosent.
    
    Som tilleggsinformasjon for å verifisere at løkken kjører riktig får vi et plot som 
    inneholder objektets bevegelse med x,y koordinater, samt hvor sluring inntreffer.
    Vi printer derfor også hvilken posisjon objektet forlot banen i i koordinater med 
    3 desimaler
    """
    
    g = 9.81 #definerer tyngdeakselerasjonen i m/s^2

    values = objects.get(object) #lagrer informasjon om et gitt objekt i en variabel
    c = values[0] #henter ut c fra dict
    m = values[1] #henter ut m fra dict
    r = sum(values[2]) #henter ut r fra dict
    v = 0 #setter startfarten lik 0
    y_0 = (R+r)*np.cos(theta_prior) #definerer startverdi for y, tar hensyn til startvinkel
    x_0 = 0 #definerer startverdi for x
    
    mu_static = 0.61 #definerer den statiske friksjonskoeffisienten
    mu_kinetic = 0.47 #definerer den kinetiske friksjonskoeffisienten
    
    energy_prior = m*g*np.cos(theta_prior)*(R+r) #startenergi med høyde justert for startvinkel

    pos_arr_y = np.zeros(2*int(1/delta_t)) 
    pos_arr_x = np.zeros(2*int(1/delta_t)) 
    #definerer to zero-arrays, et for x-posisjon og et for y-posisjon, legger også til en god
    #sikkerhetsmargin slik at løkken ikke vil bryte av at indexen ikke strekker til. 
    #Forhåndsdefinerer arrayene for at løkken skal kjøre raskt og bruke mindre minne.
    
    i = 0 #setter i for å telle antall iterasjoner
    
    f_s = 1 #definerer en f_s som er mindre enn f_s_max så løkken kjører
    f_s_max = 2 #definierer en f_s_max som er større enn f_s så løkken kjører
    
    N = 1 #definierer en vilkårlig N som er større enn 0 så løkken starter
    
    #løkke for ren rulling
    while N > 0 and f_s <= f_s_max: #løkke som kjører mens normalkraft er større enn 0, og friksjonen er mindre enn maksfriksjon
      
        alpha = ((g * np.sin(theta_prior)) / (1 + c)) / (r + R) #regner ut vinkelakselerasjon
        omega_post = omega_prior + alpha * delta_t #bruker Eulers metode for å finne vinkelfarten
        theta_post = theta_prior + omega_prior * delta_t #bruker Eulers metode for å finne vinkelen

        x = x_0 + (((R+r) * omega_prior * np.cos(theta_post) ) + v * np.cos(theta_prior)) * delta_t / 2 #regner ut x-posisjonen ved hjelp av S=1/2(v-v_0)*t, bruker v fra forrige iterasjon som v_0, og regner ut nye v direkte
        pos_arr_x[i] = x #legger til x-posisjonen i x-arrayet
        x_0 = x #oppdaterer x posisjonen
        
        y = y_0 - ((R+r) * omega_prior * np.sin(theta_post) + v * np.sin(theta_prior)) * delta_t / 2 #regner ut y-posisjonen ved hjelp av S=1/2(v-v_0)*t, bruker v fra forrige iterasjon som v_0, og regner ut nye v direkte
        pos_arr_y[i] = y #legger til x-posisjonen i x-arrayet
        y_0 = y #oppdaterer y-posisjonen
        
        v = (r + R) * omega_prior #regner ut farten
        N = m * g * np.cos(theta_prior) - (m * (v ** 2)) / (r + R) #regner ut normalkraften
        
        
        i += 1 #legger til en i antall iterasjoner
        omega_prior = omega_post #oppdaterer omega
        theta_prior = theta_post #oppdaterer theta
        f_s = c*m*alpha*(r+R) #regner ut ny friksjonskraft
        f_s_max = mu_static * N #bruker den nye normalkraften til å finne friksjonskraft
        
    break_index = i #definerer en break index når den første løkken bryter for å finne ut hvor sluring inntreffer
    
    f_k = f_s                 
    
    #løkke for rulling med sluring
    while N > 0: #løkke som kjører så lenge normalkraften er større enn 0
        alpha = ((g * np.sin(theta_prior)) - (f_k / m)) / (r + R) #regner ut vinkelakselerasjon
        omega_post = omega_prior + alpha * delta_t #bruker Eulers metode for å finne vinkelfarten
        theta_post = theta_prior + omega_prior * delta_t #bruker Eulers metode for å finne vinkelen

        x = x_0 + ( ( (R+r) * omega_prior * np.cos(theta_post) ) + v * np.cos(theta_prior)) * delta_t / 2 #regner ut x-posisjonen ved hjelp av S=1/2(v-v_0)*t, bruker v fra forrige iterasjon som v_0, og regner ut nye v direkte
        pos_arr_x[i] = x #legger til x-posisjonen i x-arrayet
        x_0 = x #oppdaterer x posisjonen
        
        y = y_0 - ((R+r)*omega_prior*np.sin(theta_post)+v*np.sin(theta_prior))*delta_t / 2 #regner ut y-posisjonen ved hjelp av S=1/2(v-v_0)*t, bruker v fra forrige iterasjon som v_0, og regner ut nye v direkte
        pos_arr_y[i] = y #legger til x-posisjonen i x-arrayet
        y_0 = y #oppdaterer y posisjonen
        
        v = (r + R) * omega_prior #regner ut hastighet
        N = m * g * np.cos(theta_prior) - (m * (v ** 2)) / (r + R) #normalkraften
        
        i += 1 #legger til en iterasjon i tellingen
        omega_prior = omega_post #oppdaterer omega
        theta_prior = theta_post #oppdaterer theta
        
        f_k = N*mu_kinetic #regner ut ny friksjonskraft

        
    
    pot_e =  m*g*(np.cos(theta_post)*(R+r)) #regner ut potensiell energi
    kin_e = 0.5 * m * v ** 2 #regner ut kinetisk energi
    
    if "hul sylinder" == object: #tar hensyn til at treghetsmomentet er ulikt for hul sylinder
        r_I = 0 
        for r_n in values[2]:
            r_I += r_n ** 2
        I = c * m * r_I #treghetsmomentet hvis vi har hul sylinder  

    else: #treghetsmoment på alle andre måter
        r_I = sum(values[2]) #henter ut radius
        I = c * m * r_I ** 2

    rot_e = 1/2 * I * (f_k*N*(R+r)*i*delta_t/I)**2 #regner ut rotasjonsenergi
        
    energy_post = pot_e + kin_e + rot_e #regner ut den totalte energien når ballen forlater banen
    
    plot_x = np.delete(pos_arr_x, [*range(i,2*int(1/delta_t))]) #lager nytt array uten nullene fra arrayet som ikke tas i bruk
    plot_y = np.delete(pos_arr_y, [*range(i,2*int(1/delta_t))]) #lager nytt array uten nullene fra arrayet som ikke tas i bruk
    
    plt.title(f"Posisjonen til {object}") #plotter tittel
    plt.xlabel("x-posisjon") #plotter navn på x-aksen
    plt.ylabel("y-posisjon") #plotter navn på y-aksen
    plt.plot(plot_x, plot_y, "r-") #plotter x og y arrayene
    
    plt.plot(plot_x[break_index], plot_y[break_index], 'bo') #plotter punktet der sluring inntreffer
    
    plt.grid() #legger til grid
    plt.show() #viser plottet
    
    func_values = f"""
{object} forlot banen etter {theta_post*180/np.pi} grader,
det er {theta_post} radianer. 

Løkken kjørte {i} ganger med delta_t: {delta_t}
{object} forlot banen etter {round(delta_t*i,6)} sekunder.

Energiendringeen ble {(energy_post - energy_prior)} joule, dette tilsvarer 
en endring i energi på {(energy_post - energy_prior) / energy_prior * 100}%

{object} forlot banen i posisjon ({round(x, 3)}, {round(y, 3)})

{energy_prior} --> {energy_post}

""" #lagrer en stor streng så svaret blir printet fint
    return func_values #returnerer strengen

print(angle(objects, "hul sylinder", R, omega_prior, alpha, theta_prior, delta_t)) #kjører og printer funksjonen
time_stop = timeit.default_timer() #stopper tiden

print(f"Programmet brukte {round(time_stop-time_start, 5)} sekunder") #printer tiden programmet brukte


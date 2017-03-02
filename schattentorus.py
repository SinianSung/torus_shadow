"""
@author         gregor.luedi
@date           2.3.2017
@description    calculating the shadow on a torus.
"""
import math
import numpy as np
import matplotlib.pyplot as plt


def sonne(time=0, majoraxis=1000, minoraxis=1003, tilt=500):
    """ Koordinaten der Sonne"""
    x_1 = majoraxis*math.cos(time)
    x_2 = minoraxis*math.sin(time)
    x_3 = tilt*math.cos(time)
    sonnenposition = np.array([x_1, x_2, x_3])
    return sonnenposition

def torus(theta, phi, radius=2, radius_tube=1):
    """Punkte auf dem Torus"""
    x_1 = (radius+radius_tube*math.cos(theta))*math.cos(phi)
    x_2 = (radius+radius_tube*math.cos(theta))*math.sin(phi)
    x_3 = radius_tube*math.sin(phi)
    point = np.array([x_1, x_2, x_3])
    return point

def normale(theta, phi, radius=2, radius_tube=1):
    """Normale zur Tangentialebene des Torus berechnen"""
    x_1 = radius_tube*(radius_tube*math.cos(phi)+radius)*math.cos(theta)*math.cos(phi)
    x_2 = radius_tube*(radius_tube*math.cos(phi)+radius)*math.sin(theta)*math.cos(phi)
    x_3 = radius_tube*(radius_tube*math.cos(phi)+radius)*math.sin(phi)
    normalenvektor = np.array([x_1, x_2, x_3])
    return normalenvektor

def diskriminante(sol, punkt):
    """ Berechnung der diskriminante """
    value_a = (np.dot(sol, sol))**2
    value_b = 4*np.dot(sol, sol)*np.dot(sol, punkt)
    value_c = 4*np.dot(sol, punkt)**2+np.dot(sol, sol)*(np.dot(punkt, punkt)+2**2-1**2)
    value_d = 4*np.dot(sol, punkt)*(np.dot(punkt, punkt)+2**2-1**2)

    return 18*value_a*value_b*value_c*value_d-4*value_b*value_b*value_b*value_d- \
        4*value_a*value_c*value_c*value_c-27*value_a*value_a*value_d*value_d

AUFLOESUNG = 50
ANZAHLSONNENSTAENDE = 12
MAJOR = 1030
MINOR = 1000
TILT = 6000

for t in range(ANZAHLSONNENSTAENDE):
    tau = 2*math.pi/ANZAHLSONNENSTAENDE*t
    sun = sonne(tau, majoraxis=MAJOR, minoraxis=MINOR, tilt=TILT)
    winkel = round(math.atan(TILT/MAJOR), 3)
    fig, ax = plt.subplots()
    plt.axis([0, 2*math.pi, 0, 2*math.pi])
    plt.title(r'Sonnenstand:  $\tau =$ {}, $\phi$: {}'.format(round(tau, 3), winkel))
    plt.xlabel(r'$p \in [0,2\pi[$')
    plt.ylabel(r'$q \in [0,2\pi[$')
    for i in range(AUFLOESUNG):
        for j in range(AUFLOESUNG):
            u = (2*math.pi)/AUFLOESUNG*i
            v = (2*math.pi)/AUFLOESUNG*j
            # Test ob der Punkt auf der sonnenabgewandten Seite liegt
            point_torus = torus(u, v)
            s = np.dot(sun, normale(u, v))
            if s < 0:
                circle = plt.Circle((u, v), 0.005, color='black')
            else:
                #Test ob der Punkt auf der entfernten Seite des Torus liegt
                temp = np.dot(sun, point_torus)
                if temp > 0:
                    circle = plt.Circle((u, v), 0.005, color='white')
                else:
                    Delta = diskriminante(sun, point_torus)
                    if Delta < 0:
                        circle = plt.Circle((u, v), 0.005, color='blue')
                    else:
                        circle = plt.Circle((u, v), 0.005, color='green')
            ax.add_artist(circle)
    fig.savefig("45/plot_{}_{}.png".format(round(tau, 3), winkel))
    plt.close()
    print('done {}'.format(t))

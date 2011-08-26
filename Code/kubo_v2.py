import sys
import os
from numpy import *

'''Program: Temp_Plot
This is a new program that I will be using to test my code. It will calculate the hall conductivity, for
the case of fermions at zero temperature. It is based off of the code temp_ptr_kubo, which is based off of
temp_ptr.py. The major change to this program, is that instead of interating over temperature, I will be
iterating over changing the fermi level. I think that this means that I will just be injecting more and more
electrons into the system, so I am not sure how the chemical potential thing will work. 

I will be able to take out all of the bisection methods however, as at zero temperature, I know exactly how the energy levels 
will be filled. I can use the simplified version of the kubo formula, as the fermi filling factor will look
like a step function.
'''



'''JOBS: 1) Get rid of degenerate energy eigenvalues- this will throw off the
probability calculations.
         2) Find the bose condensate temperature
         3) Iterate through the temperatures BEC-> END_TEMP and find the 
            chemical potential that fixes the number of particles at that
            specific temperature to a user defined value
         4) Calculate the bose occupation factor using that chemical potential
            and the temperature
         5) Output this as a file
'''


'''FUNCTION LIST'''


def n_now(work, Mu, t):
    N = 0.0
    for energy in work:
        # print t, energy, Mu
        # print exp(1)
        # energy-Mu
        # exp(1/t*(energy-Mu))
        #print "f", 1/exp(1/t*(energy-Mu))
        if energy <= t:
            n_x = 1.
        elif energy > t:
            n_x = 0.

        N = N + n_x
    #  print n_x, N
    return N
    
def bisection_t(UP, BOTTOM, n_now, energy_a, Num_Part, Mu):
    '''The purpose of this function is to calculate the bose_temperature. Which is where the temperature iteration will
    start. Mu is set to the lowest energy level'''
    #Num_part is the number of particles
    t1 = BOTTOM
    # print t1
    t2 = UP
    c = (t1+t2)/2.
    # print 'c', c
    tolerance = .00001
    for i in range(0,LARGE_NUMBER):
        
        fc = n_now(energy_a, Mu, c)-Num_Part
        f1 = n_now(energy_a, Mu, t1) -Num_Part
        f2 = n_now(energy_a, Mu, t2) - Num_Part
        #        print f1, fc, f2
        if (fc*f2 < 0):
            t1 = c
            
        elif (fc*f1 < 0):
            t2 = c
            
        else:
            print 'bisect error'
            sys.exit()
            # print 'f1',f1,'f2',f2,'fc',fc
        c = (t1+t2)/2.
            # print 'The BEC TEMP is:', t1, 'N = ', n_now(energy_a, Mu, t1),'c',c, Mu,
            
        print n_now(energy_a, Mu, t1)-Num_Part, tolerance
        if (abs(f1) < tolerance):
            print i
            break

        if i == LARGE_NUMBER:
            print "BISECTION FAILURE"
    return(t1)
        

def bisection_Mu(UP, BOTTOM, n_now, energy_a, Num_Part, temp):
    '''This is to find the right chemical potential at a specific temperature that ensures that as the temperature is
    increased, the number of particles remains conserved'''
    mu1 = BOTTOM
    # print t1
    mu2 = UP
    c = (mu1+mu2)/2
    # print 'c', c
    tolerance = .0001
    
    for i in range(0,LARGE_NUMBER):
        fc = n_now(energy_a, c, temp)-Num_Part
        f1 = n_now(energy_a, mu1, temp) -Num_Part
        f2 = n_now(energy_a, mu2, temp) - Num_Part
        print fc
        if (abs(f1) < tolerance):
            return mu1
        
        if (abs(f2) < tolerance):
            return mu2
        
        if (abs(fc) < tolerance):
            return c
        # print f1, fc, f2
        if (fc*f2 < 0):
            mu1 = c
            
        elif (fc*f1 < 0):
            mu2 = c
            
        else:
            print 'bisect error'
            sys.exit()
            #print 'f1',f1,'f2',f2,'fc',fc
        c = (mu1+mu2)/2
       #print 'The MU is:', mu1, 'N = ', n_now(energy_a, mu1, temp),'c',c, temp,

       #print n_now(energy_a, mu1, temp)-Num_Part, tolerance
    print "BISECTION FAILURE"
    sys.exit()
      
      

def readin_files(eigenName, jName):

    #The first task is to read in all the data and preserve the block structure of the matrices
    #READ-IN CURRENT OPERATOR:
    tmp1 = loadtxt(jName+'_yRE.dat', unpack=True)
    tmp2=loadtxt(jName+'_yIM.dat', unpack=True)

    tmp3 = loadtxt(jName+'_xRE.dat',unpack=True )
    tmp4=loadtxt(jName+'_xIM.dat',unpack=True )

    tpo1 = array([tmp3[0],tmp3[1]])
    jmatrix_x = vstack( (tpo1,tmp3[2::]+multiply(1j,tmp4[2::]) ))
    print jmatrix_x[1]


    tpo2 = array([tmp2[0],tmp2[1] ])
    jmatrix_y = vstack( (tpo2,tmp1[2::]+multiply(1j,tmp2[2::])) )
    jmatrix_x = jmatrix_x.transpose()
    jmatrix_y = jmatrix_y.transpose()
    #jmatrix will contain the current matrix, seperated into blocks 
    #data should be [ky][kx][ENERGY]{[E][I][G][N][E][V][E][C][T][O][R]}
    #jmatrix should be [ky,kx, m,a,t,r,i,x}
    #END READIN CURRENT OPERATOR######################

    #READ-IN HAMILTONIAN AND ENERGY#####################
    data_real = loadtxt(eigenName+'_real.dat', unpack=True)
    data_imag = loadtxt(eigenName+'_imag.dat', unpack=True)

    data  =vstack( (data_real[:4], data_real[4:]+multiply(1j,data_imag[4:])) )
    #END READ_IN HAMILTONIAN#####################
    return data, jmatrix_x, jmatrix_y


def kubo_cal(alphaVec,betaVec,alphaEnergy, betaEnergy, j_x, j_y):

    kubotop = vdot(alphaVec,dot(j_x,betaVec))*vdot(betaVec,dot(j_y,alphaVec))-vdot(alphaVec,dot(j_y, betaVec))*vdot(betaVec,dot(j_x,alphaVec))
    kubobottom = (alphaEnergy-betaEnergy+1j*1e-15)*(alphaEnergy-betaEnergy-1j*1e-15)
    return kubotop/kubobottom

            
def kubo_subband(kxList, kyList, jMatrix_x, jMatrix_y, eigMatrix ):

    #ARGS: kx is a set of kx values
    #      ky is a set of ky values
#          jMatrix_ are Q*(Num of k values)**2 matrices
#          eigMatrix is a matrix of vectors, each row contains the kx,ky, the eigenenergy and the eigenvector
#   Initialize the kubo value
    kubo=0.0
    for kx in kxList:
        for ky in kyList:#now we are in a regime where we have isolated our block. Need to test
#           Define the block as those rows which share the same kx,ky. 
            bloc =   [j[3::] for i,j in enumerate(data_trans) if data_trans[i][0] == kx and data_trans[i][1]==ky]
            bloc_jx = array([j[2::] for i,j in enumerate(jmatrix_x) if jmatrix_x[i][0] == kx and jmatrix_x[i][1] == ky])
            bloc_jy = array([j[2::] for i,j in enumerate(jmatrix_y) if jmatrix_y[i][0] == kx and jmatrix_y[i][1] == ky])

            #Now, calculate the number of points in K-space (assume Lx == Ly)
            L = float(len(jmatrix_y))
            kubo = kubo+kubo_single_energ(bloc, bloc_jx,bloc_jy,L)
    return kubo


def kubo_single_energ(Hblock, bloc_jx, bloc_jy):
    '''This function will calculate the kubo contribution from a single energy level of a subband. Ie. given a kx and a ky. 
    it is called by a function that will calculate the kubo contribution from an entire subband'''
    for alpha in Hbloc:
        for beta in Hbloc:
#           Ensure duplicates aren't called
            if all(alpha==beta):                   
                continue
                    #print bloc_jx.size, 
                    #print (len(alpha[1::])*len(beta[1::]) == bloc_jx.size)
                    #print (len(alpha[1::])*len(beta[1::]) ==bloc_jy.size), b, ky, kx
                    #print occup[alpha[0]]
                    #print occup[alpha[0]]
                    #prefac =-1j*occup[alpha[0]][-1][1]

#           test to see if the energy level is occupied on the fly ( t is the fermi energy, NOT temperature)
            if alpha[0] <= t:
                prefac = -1.0j
            else:
                prefac = 0.0
            alphaVec = alpha[1:]
            betaVec = beta[1:]
            alphaEn = alpha[0]
            betaEn = beta[0]

            
            kubo = kubo+kubo_cal(alphaVec, betaVec, alphaEn, betaEn, bloc_jx, bloc_jy)
    return kubo








#COMMAND LINE ARGS
if __name__ == '__main__':
    END_FERMI = float(sys.argv[1])
    FLUX = sys.argv[2]
    #INNAME =sys.argv[2]
    SIZEOFLATTICE = sys.argv[3]
    INTERACT = sys.argv[4]
    #INNAME2 = sys.argv[4]
    Num_Part = sys.argv[5]
    LARGE_NUMBER=1000000

    #INITIALIZE LIST VALUES##################


    NAMEBASE = 'N'+str(Num_Part)+'_Q'+str(FLUX)+'_G'+str(INTERACT)+'_Sz'+str(SIZEOFLATTICE)
    INNAME2 = 'Res/current_'+NAMEBASE
    
    INNAME = 'Res/eig_'+NAMEBASE
    Num_Part = float(Num_Part)
    print Num_Part, SIZEOFLATTICE, INNAME, INNAME2
    #
    data, jmatrix_x, jmatrix_y = readin_files(INNAME, INNAME2)
    SIZEOFLATTICE = float(SIZEOFLATTICE)
    #INIT. OCCUPTATION DICTIONARY
    occup={}
    for energy in data[3]:
        occup[energy] = []

    #INIT. KUBO LIST [t,kubo]
    kubo_list = [ ]

    #find the range of fermi energy, starting with the smallest energy level.
    beginFermi = data[3].min() 
    #Define the iteration based on the number of data points wanted
    itr = abs((END_FERMI-beginFermi)/40.0)
    t = beginFermi

    #END INITIALIZE
    contrib_file = open('./contrib_fermi.dat' ,'w')
    #iterate over the range of fermi energies. This is done by injecting more electrons into the system. N is not conserved.
    contrib_list = []
    #t_list = sort(data[3])
    while t <END_FERMI:

        #Caluculate how many electrons are present
        N = n_now(data[3], 1, t)

        data_trans = array(data.transpose())
        #First, find a list of all the different kx, ky values
        kset  = [list(set(data[1][:])), list(set(data[2][:]))]
        
        ksetj  = [list(set(jmatrix_x.transpose()[0][:])), list(set(jmatrix_x.transpose()[1][:]))]
        if set(kset[0]) != set(ksetj[0]) or set(kset[1]) != set(ksetj[1]) :
            print 'problem with k values. Check precision between two data files.'
            sys.exit()

        #The set will filter out duplicate entries, so we are left with UNSORTED, but unique k values
        kubo=0.0
        b=0
        contrib1 = 0 
        for kx in kset[0]:
            for ky in kset[1]:#now we are in a regime where we have isolated our block. Need to test
                bloc =   [j[3::] for i,j in enumerate(data_trans) if data_trans[i][1] == kx and data_trans[i][2]==ky]
                bloc_jx = array([j[2::] for i,j in enumerate(jmatrix_x) if jmatrix_x[i][0] == kx and jmatrix_x[i][1] == ky])
                bloc_jy = array([j[2::] for i,j in enumerate(jmatrix_y) if jmatrix_y[i][0] == kx and jmatrix_y[i][1] == ky])
               # print b
                #print bloc 
                trigger = 1
                contrib = 0
                for alpha in bloc:
                    for beta in bloc:
                        if all(alpha==beta):                   
                            continue
                                        #print (len(alpha)==len(beta))
                                #print bloc_jx.size, 
                                #print (len(alpha[1::])*len(beta[1::]) == bloc_jx.size)
                                #print (len(alpha[1::])*len(beta[1::]) ==bloc_jy.size), b, ky, kx
                                #print occup[alpha[0]]
                                #print occup[alpha[0]]
                                #prefac =-1j*occup[alpha[0]][-1][1]
                        if alpha[0] <= t:
                            prefac = -1.0j
                        else:
                            prefac = 0.0
                                   # assert occup[alpha[0]][-1][0] == t
                                #print alpha[1::]
                                #print bloc_jy
                                #print bloc_jx
                                #print beta[1::]
                        kubotop=(vdot(alpha[1::], dot(bloc_jx,beta[1::]))*vdot(beta[1::], dot(bloc_jy,alpha[1::]))-vdot(alpha[1::], dot(bloc_jy,beta[1::]))*vdot(beta[1::], dot(bloc_jx,alpha[1::])))
                        kubobottom=((alpha[0]-beta[0]+complex(0,1E-14))*(alpha[0]-beta[0]-complex(0,1E-14)))
                                #print prefac, kubotop, kubobottom
                               # print (dot(alpha[1::], dot(bloc_jx,beta[1::]))*dot(beta[1::], dot(bloc_jy,alpha[1::])))
                                #print kubo
                        kubo =kubo+prefac*kubotop/kubobottom/float(SIZEOFLATTICE**2)/float(FLUX)*2*pi
                        #print kubo, prefac*kubotop/kubobottom/float(SIZEOFLATTICE**2)
                        if round(prefac*kubotop/kubobottom/float(SIZEOFLATTICE**2), 4) != 0:
                            contrib = contrib +1
                            contrib1 = contrib1+1
                            contrib_list.append((t, round(prefac*kubotop/kubobottom,5 )))
                            out_line = (kubo,alpha[0], beta[0], contrib1, prefac,kubotop/kubobottom/float(SIZEOFLATTICE**2), t)
        #                    print str(out_line)
                            out_line = str(out_line)+'\n'
                            contrib_file.write(out_line)

                        #    print kubo,alpha[0], beta[0], contrib1, round(prefac*kubotop/kubobottom,5)

                                    #if int(prefac.imag) != 0:
                                     #   print prefac, prefac*kubotop/kubobottom, alpha[0], ky, kx
                                     #   if trigger != 0:
                                     #       print 'Jx_MATRIX'
                                     #       print bloc_jx
                                     #       print 'Jy_matrix'
                                     #       print bloc_jy
                                     #       print 'eig'
                                     #       print bloc
                                     #       trigger = 0




        #        print kubo, kx, ky, 'Kubo test'
                b= b+1
        kubo_list.append( (t,kubo) )
        print kubo, t, contrib1, N, SIZEOFLATTICE
        t = t+itr

    #print contrib_list[:101]
    dirname= 'PYFERMITemp'+NAMEBASE
    os.system('mkdir %s' %(dirname))
    os.system('cd ./%s' %dirname)
    os.system('pwd')
    i = 0
    out_list = occup.keys()
    outputlist= sort(out_list)

    for energy in outputlist:
        temp_a = array(occup[energy])
        name = 'F3_'+str(energy)+'target.dat'
        savetxt('%s' %dirname+'/'+name, temp_a)
        i = i+1

    savetxt('%s' %dirname+'/kubo_fermi.dat', array(kubo_list))

    contrib_file.close()
    print kubo_list

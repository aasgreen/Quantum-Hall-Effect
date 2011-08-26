import sys
import os
from numpy import *
import math 
'''Program: Temp_Plot
This is a new program that I will be using to test my code. It will calculate the hall conductivity, for
the case of bosons for the case of finite temperature. The major change to this program, is that instead of interating over temperature.

'''



'''JOBS: 
         1) Find the bose condensate temperature
         2) Iterate through the temperatures BEC-> END_TEMP and find the 
            chemical potential that fixes the number of particles at that
            specific temperature to a user defined value
         3) Calculate the bose occupation factor using that chemical potential
            and the temperature
         4) Calculate the kubo contribution and sum up over all contributions
         5) Output this as a file
'''


'''FUNCTION LIST'''


def n_now(work, Mu, t):
    N = 0.0
    for energy in work:
        #Making energy a float makes sure to cast it from a complex
        #The reason for this is because a complex will overflow too quickly 
        #So the exponential will overflow.
        energy = float(energy)
        #BOSON PREFAC##########
        if math.isinf(exp((energy-Mu)/float(t) )):
            n_x = 0.
        else:
            n_x = 1.0/(exp( (energy-Mu)/float(t))-1.0)
        #print "f", 1/exp(1/t*(energy-Mu))
#END BOSON PREFAC
        N = N + n_x
#FERMI PREFAC#####################
        #if energy <= t:
        #    n_x = 1.
        #elif energy > t:
        #    n_x = 0.
#END FERMI PREFAC

    return N
    
def bisection_t(UP, BOTTOM, n_now, energy_a, Num_Part, Mu):
    '''The purpose of this function is to calculate the bose_temperature. Which is where the temperature iteration will
    start. Mu is set to the lowest energy level. The bose temperature will be the smallest temperature where the particle number is conserved.'''
    #Num_part is the number of particles
    t1 = BOTTOM
    # print t1
    t2 = UP
    c = (t1+t2)/2.
    # print 'c', c

    #The tolerance has to be low for this to converge.
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
            
            print 'f1',f1,'f2',f2,'fc',fc
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
        #print fc
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

        #N = n_now(energy_a, mu1, temp)
        #print N,N-Num_Part, tolerance, mu1, mu2
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
    data_real = loadtxt(eigenName+'real.dat', unpack=True)
    data_imag = loadtxt(eigenName+'imag.dat', unpack=True)

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


def kubo_single_alpha(Hblock, bloc_jx, bloc_jy,alpha):
    '''This function will calculate the kubo contribution from a single energy level of a subband. Ie. given a kx and a ky. 
    it is called by a function that will calculate the kubo contribution from an entire subband'''
    kubo = 0.
    for beta in Hblock:
#       Ensure duplicates aren't called
        if all(alpha==beta):                   
            continue
            #print bloc_jx.size, 
            #print (len(alpha[1::])*len(beta[1::]) == bloc_jx.size)
            #print (len(alpha[1::])*len(beta[1::]) ==bloc_jy.size), b, ky, kx
            #print occup[alpha[0]]
            #print occup[alpha[0]]
            #prefac =-1j*occup[alpha[0]][-1][1]

#       test to see if the energy level is occupied on the fly ( t is the fermi energy, NOT temperature)
        alphaVec = alpha[1:]
        betaVec = beta[1:]
        alphaEn = alpha[0]
        betaEn = beta[0]

        
        kubo = kubo+kubo_cal(alphaVec, betaVec, alphaEn, betaEn, bloc_jx, bloc_jy)
        #print kubo_cal(alphaVec, betaVec, alphaEn, betaEn, bloc_jx, bloc_jy), alphaEn, betaEn, 'kubo_single_alpha'
    return kubo


def occupation(record,energyVec,t,mu,n_now):
    '''Calculate the occupation and probability of a specific energy level'''
    '''ARGS: record is the ongoing record of the ocuppation as a function of t. It has
        type dictionary, and must be initialized before this function is called!
        energyVec is a vector containing all the energy states possible
        t is the temperature
        mu is the chemical potential
        n_now is a function that calculates how many particles there are.
        '''
    testN = 0.
    testP = 0.
    for energy in energyVec:
        energy = energy.astype(float)
        if math.isinf(exp( (energy-float(mu))/float(t) )):
            ocup = 0.
            print 'inf', exp( (energy-mu)/t), math.isinf( exp( (energy-mu)/t) )
        else:
            ocup = 1./(exp( (energy-float(mu))/float(t) )-1.)
         
        N = n_now(energyVec, mu, t)
        prob = ocup/N
        #Check for degeneracies
        if record[energy] != []:
            if record[energy][-1] == [t,ocup,prob]:
                record[energy][-1] = [t,record[energy][-1][1]+ocup, record[energy][-1][2]+prob]
                print record[energy][-1]
            else:
                record[energy].append([t,ocup, prob])

        else:
            record[energy].append([t,ocup,prob])
        testN = testN+ocup
        testP = testP +prob
    #print testN, N, testP
    assert( testN == N)
    
       
    return record


def basis_states(N, SIZEOFLATTICE, FLUX):
#We need to run the process basis generator, to find out how
#many ways needed to distribute N particles amongst the number of states
#available.

#Calculate the number of states avaiable to the system.
    N = int(N)
    nStates = SIZEOFLATTICE**2*FLUX
    nStates = int(nStates)
#Run the process basis_gen.sh, which takes as arguments: Number of particles, Number of states
#and will generate a file called inff1.dat, which can be read in, this will contain the enumerate
#different combos/perms of N particles in nStates.

    run_cmd = './basis_gen.sh '+str(N)+' '+str(nStates)
    os.system(run_cmd)

#Now read in the resulting file:

    bStates = loadtxt('inff1.dat')
    os.system('rm inff1.dat')
#be a good citizen
    return bStates


def main_kubo(bState,data,jmatrix_x, jmatrix_y, SIZEOFLATTICE, FLUX,mu,t,N):
    kubo = 0.
    for basis in bState: #this sums over the 'alphas' in the kubo formula
        #basis will give way to combine block eigenvectors to get the alpha
        alphaEnergy ={}
        kubo_state = {}
        #Split the basis into blocks, each Q long, knowing total length = Q*SIZEOFLATTICE**2
        basisBlocks = hsplit(basis, SIZEOFLATTICE**2)
        for startIndex,block in enumerate(basisBlocks):
            #Each block is isolated from the others, so we can find the contribution from each
            #of these.
            numberParticlesInBloc = sum(block)
            #print numberParticlesInBloc
            if numberParticlesInBloc != 0:
                #Then this block will contribute to kubo

                sizeOfBand = len(data[numberParticlesInBloc][0])/SIZEOFLATTICE**2
                Hblock = data[numberParticlesInBloc][3:].T[startIndex*sizeOfBand:(startIndex+1.)*sizeOfBand]
                
                jxBlock = jmatrix_x[numberParticlesInBloc].T[2:].T[startIndex*sizeOfBand:(startIndex+1.)*sizeOfBand]
                jyBlock = jmatrix_y[numberParticlesInBloc].T[2:].T[startIndex*sizeOfBand:(startIndex+1.)*sizeOfBand]
                
                #Find the unique basis state vector associated with the energy
                alphaEnergy[startIndex] = 0.
                for blockIndex,stateValue in enumerate(block):
                    alphaEnergy[startIndex] = alphaEnergy[startIndex]+stateValue*float(data[1][3][startIndex*FLUX+blockIndex])
             #       print data[1][3][startIndex*FLUX+blockIndex], 'single particle energy', stateValue
            #    print alphaEnergy[startIndex], 'total energy of block'

                alphaVec = [i for i in Hblock if round(i[0].astype(float)-alphaEnergy[startIndex],12) == 0] 

           #     print alphaVec, 'alphaVec'
           #     print Hblock[-1].astype(float), 'Hblock'
           #     print alphaEnergy[startIndex],Hblock[-1][0].astype(float)
           #     print alphaEnergy[startIndex] == Hblock[-1][0].astype(float)
          #      print alphaEnergy[startIndex]- Hblock[-1][0].astype(float)
                assert size(alphaVec) == size(Hblock[0])
                #if this fails, then we have degeneracies and this is no longer a
                alphaVec = array(alphaVec).squeeze()
                #way to uniquely identify the eigenstates

#                print jyBlock.T[1:].T, 'jy'
                kubo_state[startIndex] = kubo_single_alpha(Hblock, jxBlock, jyBlock, alphaVec)
         #       print kubo_state[startIndex],'kubo_state[startIndex]'
       #end for basisblocks
       #Each basisblock gave a kubo contribution.
       #Should now have a list of kubo states. These will be addeded together and then
       #the entire term will be weighted by the total energy
        total_energy = sum([alphaEnergy[i] for i in alphaEnergy.keys()])
        #print total_energy, 'total energy'
        total_particle = sum(basis)
        prefac = -1.j/(exp( (total_energy-mu)/t )-1)/N
        #print prefac, 'prefac', exp( (total_energy -mu)/t), mu, t, exp( (total_energy-mu*total_particle)/t)-1, total_particle
        #This 'big alpha' is the summation of specific combination of eigenstates of blocks. Now
        #iterate through them and add up the contributions from each block to the kubo. Then weight
        #by the total energy of the alpha state.
        for contribution in kubo_state.keys():
                kubo = kubo + prefac*kubo_state[contribution]/SIZEOFLATTICE**2/float(FLUX)*(2*pi)
    #    print kubo, 'this is the kubo', kubo/prefac*SIZEOFLATTICE**2
           #Iterate onto the next alpha value...
    return kubo

#COMMAND LINE ARGS
if __name__ == '__main__':
    END_TEMP = float(sys.argv[1])
    FLUX = sys.argv[2]
    SIZEOFLATTICE = sys.argv[3]
    INTERACT = sys.argv[4]
    Num_Part = sys.argv[5]
    LARGE_NUMBER=10000
    #INITIALIZE LIST VALUES##################
    NAMEBASE = 'N'+str(Num_Part)+'_Q'+str(FLUX)+'_G'+str(INTERACT)+'_Sz'+str(SIZEOFLATTICE)
    INNAME2 = 'Res/current_'+NAMEBASE
    
    INNAME = 'Res/eig_'+NAMEBASE+'_'
    Num_Part = float(Num_Part)
    SIZEOFLATTICE=float(SIZEOFLATTICE)
    FLUX = float(FLUX)
    #CALCULATE BASIS STATE OF ENTIRE SYSTEM:
    bState = basis_states(Num_Part, SIZEOFLATTICE, FLUX)

    #READIN all data associated with the different particle states.
    data = {}
    jmatrix_x = {}
    jmatrix_y = {}
    for number in range(1,int(Num_Part+1)):
        numPartFileName = 'N'+str(number)+'_Q'+str(int(FLUX))+'_G'+str(INTERACT)+'_Sz'+str(int(SIZEOFLATTICE))
        INNAME2number = 'Res/current_'+numPartFileName
        
        INNAMEnumber = 'Res/eig_'+numPartFileName+'_'
 
        data[number], jmatrix_x[number], jmatrix_y[number] = readin_files(INNAMEnumber, INNAME2number)
    #print size(bState[0]), size(data[1][3])
    fullSpectrum = dot(bState, data[1][3])
    #INIT. OCCUPTATION DICTIONARY
    occup={}
    for energy in fullSpectrum:
        occup[energy] = []

    #INIT. KUBO LIST [t,kubo]
    kubo_list = [ ]
    #find the range of fermi energy, starting with the smallest energy level.
    lowestEnergy = fullSpectrum.min()
    print 'The lowest energy is', lowestEnergy
    #Define the iteration based on the number of data points wanted
    itr = END_TEMP/20
    t = lowestEnergy


    #SOLVE FOR BOSE TEMPERATURE###############
    Mu_UPBOUND = float(lowestEnergy)-.00001
    
    #So, the problem is finding the value of t that will give the
    #proper N == Num_Part at a Mu= MU_UPBOUND
    #We know that when MU == MU_UPBOUND, we are at the beginning of
    #a bose einstien condesate. The number of particles should still
    #be equal to Num_part, as any temp lower than T_BEC, the chemical
    #potential will have to be fixed at MU_UPBOUND, so the particle
    #number will have to start decresing. 

    #We can just use the bisection method, as we know that the temperature
    #will be very small. Just choose an arbitrarily large number as one 
    #bound and zero as another
    UP = 10.0
    BOTTOM = -1.0
    t_begin = bisection_t(UP, BOTTOM, n_now, fullSpectrum, Num_Part, Mu_UPBOUND)

    #Now, perturb from this, to avoid rounding errors
    t = t_begin+float(t_begin)/1000.
    #END SOLVE BOSE TEMPERATURE###########


    #To speed up the bisection, I can calculate the maximum possible mu, this will give a slightly
    #smaller area to work with.
    MAX_MU = bisection_Mu(-10E22, Mu_UPBOUND, n_now, fullSpectrum, Num_Part, END_TEMP+1)

    #Initialize the dictionary that will hold the occupation factors as a function of time.
    ocup = {}
    for energy in set(fullSpectrum):
        occup[energy] = [[]]
    #END INITIALIZE##############

    #Begin MAIN LOOP OVER TEMPERATURE##################################
    while t < END_TEMP:
        #First, calculate the mu, fixing the number of bosons present
        mu = bisection_Mu(MAX_MU, Mu_UPBOUND, n_now, fullSpectrum, Num_Part, t)
        Mu_UPBOUND = mu
       	N = n_now(fullSpectrum, mu, t)
        assert round(N,3) == float(Num_Part)
        kubo = main_kubo(bState, data, jmatrix_x, jmatrix_y, SIZEOFLATTICE, float(FLUX),mu,t,N)
        kubo_list.append( (t,kubo) )
#        occup = occupation(occup, fullSpectrum, t, mu, n_now)
        print kubo,t ,'kubo/temp'
        t = t+itr

    dirname= 'PYBOSONTemp_'+str(END_TEMP)+'_'+NAMEBASE
    os.system('mkdir %s' %(dirname))
    os.system('cd ./%s' %dirname)
    os.system('pwd')
    savetxt('%s' %dirname+'/kubo_fermi.dat', array(kubo_list))




    '''
    chaff
 bloc =   [j[3::] for i,j in enumerate(data_trans) if data_trans[i][1] == kx and data_trans[i][2]==ky]
                                bloc_jx = array([j[2::] for i,j in enumerate(jmatrix_x) if jmatrix_x[i][0] == kx and jmatrix_x[i][1] == ky])
                                bloc_jy = array([j[2::] for i,j in enumerate(jmatrix_y) if jmatrix_y[i][0] == kx and jmatrix_y[i][1] == ky])
                 for alpha in block[m]:
                     for beta in block[m]:
                         kubo=0.0
                         if all(alpha==beta):                   
                             continue
                         kubotop=(vdot(alpha[1::], dot(bloc_jx,beta[1::]))*vdot(beta[1::], dot(bloc_jy,alpha[1::]))-vdot(alpha[1::], dot(bloc_jy,beta[1::]))*vdot(beta[1::], dot(bloc_jx,alpha[1::])))
                         kubobottom=((alpha[0]-beta[0]+complex(0,1E-14))*(alpha[0]-beta[0]-complex(0,1E-14)))
                         kubo_state[m] =kubo_state[m]+prefac*kubotop/kubobottom/float(SIZEOFLATTICE**2)
                         if round(prefac*kubotop/kubobottom/float(SIZEOFLATTICE**2), 4) != 0:
                             contrib = contrib +1
                             contrib1 = contrib1+1
                             out_line = (float(kubo),alpha[0], beta[0], contrib1, prefac,kubotop/kubobottom/float(SIZEOFLATTICE**2),t)
                             print str(out_line)
                             out_line = str(out_line)+'\n'
                             contrib_file.write(out_line)
                             b= b+1
                             '''

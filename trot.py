#! /usr/bin/env python
import numpy as np
import subprocess

# See example at the bottom of this module

class Trot:
    """
    Class providing functionality to rotate a tilted-ring model as defined in a tirific.def file

    Instance variables:
        randstate (numpy random object): status of random generator (see __init__)
        tirexec (string): tirific executable

    Methods:
        __init__      Initialise random number generator
        readdef       Read a tirific deffile and return a list with parameters
        writedef      Write a deffile based on parvallist, replacing PA and INCL by pa and incl (which are assumed to be a list of floats)
        getline       Return the value of a parameter as numpy darray
        normfrompincl Convert a numpy darray of position angles and a numpy darray of inclinations into a darray of corresponding normal vectors
        rotnv         Rotate darray of vectors nv each around z-axis by alpha, then around x-axis by beta, then around z-axis by gamma
        pinclfromnv   Convert a 3xn darray on n vectors into a list of position angles and a list of inclinations
        getralph      Calculate random numbers two of which can be interpreted as being inclination and position angle of normal vectors randomly distributed across a unit sphere
    """
    
    def __init__(self, seed = None, tirexec = 'tirific'):
        """
        Initialise random number generator

        Variables:
        seed (any type) Seed of the random number generator. If not long or int, a random-generated number (see doc of numpy.random
        Init sets the pseudo random number generator (seed may be an integer to make the process deterministic, otherwise the random generator is initiated by an arbitrary number. Instance variable
        """
        np.random.seed(seed)
        self.randstate = np.random.get_state()
        self.tirexec = tirexec

    def readdef(self, deffile):
        '''
        Read a tirific deffile and return a list with parameters

        Parameters:
        deffile (string): input deffile name

        Return: readdef (list) List of pairs with parameter name (string) and value (string)

        Opens deffile and reads content into a list of pairs
        consisting of parameter name and value, then closes
        deffile. This is a primitive realisation, recognising the = -
        symbol as separator and assuming that all parameter - value
        pairs are separated by carriage returns.

        '''
        file = open(deffile)
        parvallist = []
        for line in file.readlines():
            parvallist += [line.split('=')]
            file.close()
        return parvallist

    def replaceparvallist(self, parvallist, par, val):
        """Replace parameter values of par with values val, converted to string and space as delimiter
        
        Input:
        parvallist (list)  : List of pairs consisting of parameters (string) and their values (string)
        par (string)       : Parameter name
        val (list of float): Parameter values

        Scans parameter names in parvallist until it finds name
        (ignoring whitespaces) and replaces the corresponding value
        with val, converted to a string
        """
        for line in parvallist:
            if line[0].strip() == par:
                thestring = ''
                for number in val:
                    thestring += ' %.6E' % number
                line[1] = thestring+'\n'
        return

    
    def writedef(self, deffile, parvallist):
        '''
        Write a deffile based on parvallist

        Parameters:
        deffile (str):        Output tirific .def file
        parvallist (list):    List of pairs consisting of parameters (string) and their values (string)

        Return: void

        Write every pair in parvallist to deffile, adding a = - symbol
        in-between.
        '''
        file = open(deffile, 'w')
        for line in parvallist:
            if len(line) == 0:
                pass
            elif len(line) == 1:
                file.write(line[0])
            else:
                file.write(line[0]+'='+line[1])
        file.close()
        return

    def getline(self, parvallist, parametername):

        """
        Return the value of a parameter as numpy darray

        Input:
        parvallist (list of pairs of string): List of parameter-value pairs
        parametername (string):               Name of parameter

        Scans the parameter names in parvallist until parametername is
        found and split the value of that parameter by space to return
        the resulting list converted to a numpy double darray.

        """
        for line in parvallist:
            if line[0].strip() == parametername:
                if line[0].strip() == parametername:
                    try:
                        a = line[1].split()
                    except:
                        return np.array([])
                    try:
                        b = np.array(a, dtype='double')
                    except:
                        return np.array([])
        return b

    def normvfrompincl(self, pa, incl):
        """
        Convert a numpy darray of position angles and a numpy darray of inclinations into a darray of corresponding normal vectors

        Input:
        pa (double darray):   Darray of position angles in radians
        incl (double darray): Darray of inclinations in radians
        
        Return:
        normfrompincl: 3xn darray of n normal vectors

        Calculate the normal vectors corresponding to the n pairs pa
        and incl and return a 3xn darray of those normal vectors.

        """
        return np.stack((-np.sin(pa)*np.sin(incl), np.cos(pa)*np.sin(incl), np.cos(incl)))

    def rotnv(self, nv, alpha, beta, gamma):
        """
        Rotate darray of vectors nv each around z-axis by alpha, then around x-axis by beta, then around z-axis by gamma
        
        Input:
        nv (float darray): 3xn array of vectors
        alpha (float):     angle to rotate around z-axis in radians
        beta (float):      angle to rotate around x-axis in radians
        gamma (float):     angle to rotate around z-axis in radians

        Return: 3xn-darray of n rotated vectors
        
        Rotate each of the n vectors in the input array first by alpha
        around z-axis, then around x-axis by beta, then around z-axis
        by gamma, return the result.

        """
        alpar = [[np.cos(alpha), -np.sin(alpha), 0.],[np.sin(alpha), np.cos(alpha), 0.],[0., 0., 1.]]
        betar = [[1., 0., 0.], [0., np.cos(beta), -np.sin(beta)],[0., np.sin(beta), np.cos(beta)]]
        gamar = [[np.cos(gamma), -np.sin(gamma), 0.],[np.sin(gamma), np.cos(gamma), 0.],[0., 0., 1.]]
        return np.dot(gamar,np.dot(betar,np.dot(alpar,nv)))

    def pinclfromnv(self, nv):
        """
        Convert a 3xn darray on n vectors into a darray of position angles and a darray of inclinations

        Input:
        nv (float darray): 3xn darray representing n normal vectors of rings

        Return: pinclfromnv pair consisting of a darray of position angles and a darray of inclinations, both in radians

        Takes as input a 3xn float darray of n vectors which are
        interpreted to be normal vectors of rings, each of which with
        a position angle and an inclination, which are returned each
        as a darray. The method is a bit ugly to avoid floating errors
        when dealing with small numbers.

        Notice that of course position angle and inclination pairs are
        ambiguous, e.g. incl = 350 pa = 50 describes the same normal
        vector as incl=10 pa=230. The function will always return 0 <=
        incl < 180 deg (pi) and 0 <= incl < 360 deg (2pi).

        """
        x = nv[0]
        y = nv[1]
        z = nv[2]
        pa = np.mod(np.arctan2(-x, y), 2*np.pi)
        # Instead of using incl = np.arctan2(y/np.cos(pa),z) or incl = np.arctan2(-x/np.cos(pa),z), we use both ways to calculate incl, in regions were cos and sin respectively are not close to 0. A bit ugly...

        pac = 4.*pa/np.pi
        conditions = [(0. <= pac) & (pac < 1), (1 <= pac) & (pac < 3), (3 <= pac) & (pac < 5), (5 <= pac) & (pac < 7), (7 <= pac) & (pac < 8)]
        f1 = lambda x: 1./np.cos(x)
        f2 = lambda x: -1./np.sin(x)
        produ = np.piecewise(pa, conditions, [f1, f2, f1, f2, f1])
        g1 = lambda x: 0.
        g2 = lambda x: x
        xsel = np.piecewise(x, conditions,[g1, g2, g1, g2, g1])
        ysel = np.piecewise(y, conditions,[g2, g1, g2, g1, g2])
        xy = (xsel+ysel)*produ
        incl = np.arctan2(xy,z)
        return pa,incl
        
    def getralph(self):
        '''
        Calculate random numbers two of which can be interpreted as being inclination and position angle of normal vectors randomly distributed across a unit sphere

        Input:
    
        Return: getralph list of angles

        Return three random angles alpha, beta, and gamma, alpha
        uniformly between 0 and 2 pi, beta between -pi/2 and pi/2,
        such that beta is distributed like sin beta, and gamma
        uniformly between 0 and 2 pi. These can be used to describe a
        rotation by alpha about the z-axis followed by a rotation
        about the x-axis by beta and about the z-axis by gamma. This
        way, if the normal vector is initially on the z-axis, it is
        uniformly distributed on a unit sphere after the rotations.
        The random number status will be read from the instance
        variable randstate and written to the same variable at the
        beginning and the end of the function respectively.

        '''
        np.random.set_state(self.randstate)
        alpha = np.random.uniform(0, 2*np.pi)
        beta = np.arcsin(np.random.uniform(-1,1))
        gamma = np.random.uniform(0,2*np.pi)
        self.randstate = np.random.get_state()
        return [alpha, beta, gamma]
    
    def rotatefrom(self, parvallist, paname, inclname, alpha, beta, gamma):
        """
        Take inclination with name inclname and position angle with name paname and rotate using alpha, beta, and gamma

        Input:
        parvallist (list): List of pairs consisting of parameters (string) and their values (string)
        paname (string):   Name of string representing position angles in parvallist (e.g. PA, PA_2, PA_3, ...)
        inclname (string): Name of string representing inclinations in parvallist (e.g. INCL, INCL_2, INCL_3, ...)
        alpha (float):     Angle to rotate around z-axis in radians
        beta (float):      Angle to rotate around x-axis in radians
        gamma (float):     Angle to rotate around z-axis in radians

        Output: void

        Find values corresponding to paname and inclname in parvallist
        and interpret them as inclination and position angle. Convert
        into a normal vector (0,0,1) and rotate the normal vector
        first around z-axis by alpha, then around x-axis by beta, then
        around z-axis by gamma. Then convert back to position angle
        and inclination and replace the values in parvallist by the
        new ones (in degrees).

        """
        pa = np.pi*self.getline(parvallist, paname)/180.
        incl = np.pi*self.getline(parvallist, inclname)/180.
        nv = self.normvfrompincl(pa, incl)
        newnv = self.rotnv(nv, alpha, beta, gamma)
        pa, incl = self.pinclfromnv(newnv)
        pa = 180.*pa/np.pi
        incl = 180.*incl/np.pi
        self.replaceparvallist(parvallist, paname, pa)
        self.replaceparvallist(parvallist, inclname, incl)
        return

    def rotatesingle(self, indeffile, outdeffile, outcube, alpha, beta, gamma):
        """
        Rotate a tilted-ring model as specified by indeffile and produce an output cube

        Input:
        indeffile (string): Input deffile name
        outdeffile (str):   Output tirific .def file name
        outcube (str):      Output cube name
        alpha (float):     angle to rotate around z-axis in radians
        beta (float):      angle to rotate around x-axis in radians
        gamma (float):     angle to rotate around z-axis in radians

        Return: void

        Takes indeffile as an input and rotates all inclinations and
        position angles occuring (INCL, INCL_2, ..., PA, PA_2, ...)
        first by alpha around z-axis (LOS), then by beta about x-axis
        (N), then by gamma about z-axis (LOS). Writes result into
        outddeffile and runs tirific to return outcube (blocking any
        other output). Caution: rotation is done by converting
        node-wise to normal vector, then rotating normal vector, then
        converting back to inclination and position angle. Ambiguities
        may result in discontinuities.

        """
        parvallist = self.readdef(indeffile)
        ndisks = int(self.getline(parvallist, 'NDISKS')[0])
        self.rotatefrom(parvallist, 'PA', 'INCL', alpha, beta, gamma)
        for i in range(2,ndisks+1):
            self.rotatefrom(parvallist, 'PA_%i' % i, 'INCL_%i' % i, alpha, beta, gamma)
        print 'Generating deffile %s' % outdeffile
        self.writedef(outdeffile, parvallist)
        command = self.tirexec+' LOOPS= 0 GR_CONT= ACTION= 1 LOGNAME= TABLE= BIGTABLE= TIRDEF= TIRSMO= COOLGAL= TILT= BIGTILT= INCLINO= GR_DEVICE= DEFFILE='+outdeffile+' OUTSET='+outcube
        #print command
        print 'Generating model %s' % outcube
        output = subprocess.Popen(command.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        #print output.communicate()
        return

    def rotatemulti(self, indeffile, outdefprefix, outcubeprefix, nmodels):
        """
        Rotate a tilted-ring model to random direction and produce an output cubes

        Input:
        indeffile (string) : Input tirific deffile name
        outdefprefix  (str): Prefix of output tirific .def file names
        outcubeprefix (str): Output cube name prefixes
        nmodels (int)      : Number of output models

        Return: void

        Takes indeffile as an input and generates nmodels random
        rotations by means of triplets of rotation angles (alpha,
        beta, gamma) to rotate all inclinations and position angles
        occuring (INCL, INCL_2, ..., PA, PA_2, ...)  first by alpha
        around z-axis (LOS), then by beta about x-axis (N), then by
        gamma about z-axis (LOS). Writes result into outddefprefix_i.def,
        where i is the number of the output model and runs tirific to
        return outcubeprefix_i.fits (blocking any other output). Caution:
        rotation is done by converting node-wise to normal vector,
        then rotating normal vector, then converting back to
        inclination and position angle. Ambiguities may result in
        discontinuities.

        """

        # This can certainly be done in a more elegant way
        for j in range(nmodels):
        
            # n random rotations with actual generation of cubes
            alpha, beta, gamma = self.getralph()
            parvallist = self.readdef(indeffile)
            ndisks = int(self.getline(parvallist, 'NDISKS')[0])
            self.rotatefrom(parvallist, 'PA', 'INCL', alpha, beta, gamma)
            for i in range(2,ndisks+1):
                self.rotatefrom(parvallist, 'PA_%i' % i, 'INCL_%i' % i, alpha, beta, gamma)
            outdefname = outdefprefix+'_%i.def' % j
            print 'Generating deffile %s' % outdefname
            self.writedef(outdefname, parvallist)
            outcubename = outcubeprefix+'_%i.fits' % j
            command = self.tirexec+' LOOPS= 0 GR_CONT= ACTION= 1 LOGNAME= TABLE= BIGTABLE= TIRDEF= TIRSMO= COOLGAL= TILT= BIGTILT= INCLINO= GR_DEVICE= DEFFILE='+outdefname+' OUTSET='+outcubename
            #print command
            print 'Generating model %s' % outcubename
            output = subprocess.Popen(command.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            #print output.communicate()
        return
        
#if __name__ == '__main__':
#    return

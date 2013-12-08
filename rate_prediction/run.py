_CVSVers="$Id: sram4631.py,v 1.21 2010/03/31 16:50:24 sierawbd Exp $"

import sys
from optparse import OptionParser
import base64
import cPickle
import G4Support
import math
import cPickle
import string

###############################################################################
# Extra Units
###############################################################################
nm=1e-6*mm
meter=1000*mm

# A list of semiconductor devices
deviceList=[]
mbuPatterns=[]

def PostProcess():
    LET=None
    if not mred.gun.energy_spectrum:
        LET=mred.LET(mred.gun.particle,'silicon',mred.gun.energy)
        density=mred.materials.materials_dict['silicon'].GetDensity() / (g / cm/cm/cm)
        print "Finished running %0.3g MeV %s (LET=%0.3g MeV-cm2/mg) at %0.1f az, %0.1f el" % (mred.gun.energy, mred.gun.particle.GetParticleName(), LET, options.azimuth, options.elevation)
    # Add some additional attributes to the HDF5 files. With these values,
    # postprocessing scripts can easily create reverse integral cross section
    # or event rate plots.   
    if mred.hdf5.write_output_files:
        if mred.gun.random_spatial_sampling == "directionalFlux":
            mred.hdf5.setFileAttribute('fluence_unit', mred.gun.fluence_unit)
            #mred.hdf5.setFileAttribute('beamTilt', options.beamTilt)
            #mred.hdf5.setFileAttribute('beamRoll', options.beamRoll)        
        rd=mred.detCon.GetDeviceRadius()/cm
        mred.hdf5.setFileAttribute('rd', rd)
        mred.hdf5.setFileAttribute('nIons', mred.runMgr.GetEventCount())
        if LET: mred.hdf5.setFileAttribute('LET', LET)
        mred.hdf5.setFileAttribute('nDevices', (options.nRows*options.nCols))
        #mred.hdf5.setFileAttribute('args', run_args[1:])

        if globals().has_key('_CVSVers'):
            mred.hdf5.setFileAttribute('CVStag', _CVSVers)

    if len(mbuPatterns) > 0:
        # Dump the upset patterns to a file for postprocessing
        print "Outputing %d patterns" % (len(mbuPatterns))
        fh=file(pickleFile,"w")
        cPickle.dump(mbuPatterns,fh)
        fh.close()

    if options.spectrum and options.particle == "neutron":
        flux=0
        lowFlux=0
        highFlux=0
        steps=1000
        minLogE=math.log10(mred.gun.energy_spectrum.low_energy)
        maxLogE=math.log10(mred.gun.energy_spectrum.high_energy)
        dLogE=(maxLogE-minLogE)/float(steps)
        for i in range(steps):
            E=10**(minLogE+i*dLogE)
            dE=10**(minLogE+(i+1)*dLogE) - 10**(minLogE+i*dLogE)
            if E >= mred.gun.energy_spectrum.low_energy:
                flux += mred.gun.energy_spectrum.spectrum(E) * dE
            if E >= 1 and E > mred.gun.energy_spectrum.low_energy:
                lowFlux += mred.gun.energy_spectrum.spectrum(E) * dE
            if E >= 10 and E > mred.gun.energy_spectrum.low_energy:
                highFlux += mred.gun.energy_spectrum.spectrum(E) * dE
        if mred.hdf5.write_output_files:
            mred.hdf5.setFileAttribute('total_flux', flux)
            mred.hdf5.setFileAttribute('flux_above_tenMeV', highFlux)
        print "Total flux %0.3g n/cm2/h" % ( flux * 3600 )
        print "Flux > 1.0 MeV %0.3g n/cm2/h" % ( lowFlux * 3600)
        print "Flux > 10.0 MeV %0.3g n/cm2/h" % ( highFlux *3600 )
        fluence = mred.gun.fluence_unit * mred.runMgr.GetEventCount()
        print "Fluence = %0.3g n/cm2" % ( fluence )
        print "Valid = %g %g" % ( validEvents_y, validEvents_y2)
        yscaler = (1e9 * 60 * 60 * 2**20)/fluence
        print "FIT = %g +- %g" % ( validEvents_y * yscaler, math.sqrt(validEvents_y2) * yscaler )

###############################################################################
# Parsing Options
###############################################################################
parser=OptionParser(usage="myScript [options] [script files]")
parser.disable_interspersed_args()
# Options for run
parser.add_option("", "--suv", action='store_true', dest="enableSUV", default=False,
                  help="Use the GEANT4 OpenGL viewer")
parser.add_option("", "--dx", action='store_true', dest="enableDX", default=False,
                  help="Use the OpenDX viewer for event-by-event viewing")
parser.add_option("", "--dx-export", action='store_true', dest="exportDX", default=False,
                  help="Use the OpenDX viewer to visulize the target")
parser.add_option("", "--exportTracks", action='store_true', dest="exportTracks", default=False,
                  help="Write the tracks of valid events to HDF5")
parser.add_option("", "--nIons", action='store', dest="nIons", type="int", default=100,
                  help="The number of particles to run")
parser.add_option("", "--particle", action='store', dest="particle", type="str", default="proton",
                  help="Particle species")
parser.add_option("", "--beamE", action='store', dest="beamE", type="float", default=100.0,
                  help="Beam energy (MeV)")
parser.add_option("", "--beamP", action='store', dest="beamP", type="float", default=None,
                  help="Beam momentum (MeV/c)")
parser.add_option("", "--beamZ", action='store', dest="beamZ", type="int", default=1,
                  help="Beam atomic number")
parser.add_option("", "--beamA", action='store', dest="beamA", type="int", default=1,
                  help="Beam atomic weight")
parser.add_option("", "--spectrum", action='store', dest="spectrum", type="str", default='/home/sierawbd/environments/PARMA_sea_level_nyc/PARMA_mu_minus.dat',
                  help="Radiation energy spectrum")
parser.add_option("", "--environment", action='store', dest="environment", type="str", default=None,
                  help="Radiation environment")
parser.add_option("", "--event", action='store', dest="random_seeds", type="str", default=None,
                  help="Replay the event given by the random seeds")
# Options for physics
parser.add_option("", "--no-recoils", action='store_true', dest="disableRecoils", default=False,
                  help="Enable/disable energy deposition from scattering events")
parser.add_option("", "--no-nuclear", action='store_true', dest="disableNuclear", default=False,
                  help="Enable/disable energy deposition from scattering events")
parser.add_option("", "--sigmaBias", action='store', dest="sigmaBias", type="float", default=200,
                  help="Hadronic cross section bias factor")
parser.add_option("", "--rangeCuts", action='store', dest="rangeCuts", type="float", default=1.0,
                  help="Range cuts for electrons (um)")
# Options for batch submission
parser.add_option("", "--runName", action='store', dest="runName", type="str", default=None,
                  help="Simulation run name")
# Options for target structure
parser.add_option("", "--rows", action='store', dest="nRows", type="int", default=10,
                  help="The number of device rows")
parser.add_option("", "--cols", action='store', dest="nCols", type="int", default=10,
                  help="The number of device columns")
parser.add_option("", "--tech", action='store', dest="tech", type="float", default=65.0,
                  help="The scaling factor for the sensitive volume model")
parser.add_option("", "--process", action='store', dest="process", type="str", default="bulk",
                  help="The process for the sensitive volume model")
parser.add_option("", "--min_energy", action='store', dest="min_energy", type="float", default=0,
                  help="The minimum energy deposited in a device")

(options, args) = parser.parse_args(run_args[1:])
print "BATCH_VARS", dir(batch_vars)
print "RUN_ARGS", run_args

###############################################################################
# Configure Physics
###############################################################################
PyG4Core.ExecuteCommand('/control/verbose 2')
#PyG4Core.ExecuteCommand('/control/saveHistory')
PyG4Core.ExecuteCommand('/run/verbose 2')
if mred.physics.list_name == 'PhysicsList':
    # --------------------------------------------------------------------------
    # 'EmStandard' is the c++ default EM model. Select from this list if another
    # module is desired. Selecting 'EmStandard' here is unnecessary. It is listed
    # only for completeness.
    mred.physics.addModule('EmStandard')			# The c++ default. Same as G4 EM option 3.
    # mred.physics.addModule('EmPenelopeQED')			# This is recommended for most applications.
    # mred.physics.addModule('EmStandardScreened')	# Alternative general module.
	
    # mred.physics.addModule('EmLowEnergyQED')		# The VU Livermore low energy physics builder.
    # mred.physics.addModule('G4EmStandard_opt1')	# G4 builder.
    # mred.physics.addModule('G4EmStandard_opt2')	# G4 builder.
    # mred.physics.addModule('G4EmStandard_opt3')	# G4 builder. Same as 'EmStandard'.
    # mred.physics.addModule('EmLowEnergy')			# Geant4 Livermore builder.
    # mred.physics.addModule('G4EmLivermorePolarizedPhysics') # G4 builder.
    # mred.physics.addModule('EmPenelope')			# The Geant4 Penelope builder.
    #
    # --------------------------------------------------------------------------
    # These are for elementary particles. Normally use them all.
    if not options.disableNuclear:
        mred.physics.addModule('HadronElastic')
        mred.physics.addModule('HadronInelastic')
        mred.physics.addModule('PiKInelastic')
    # --------------------------------------------------------------------------
    # Neutron & proton hadronic physics (elastic and inelastic). This is the
    # standard physics module for protons and neutrons. It contains a number
    # of switches to control behavior, but the standard configuration
    # for protons and neutrons is to use CEM03 up to 5.5 GeV/u and QGSC
    # from 5 GeV/u to 100 TeV/u. The two overlap in the range 5-5.5 GeV/u.
    # By default high precision neutron physics is used below 20 MeV/u. The 
    # switches summarized below with their default values can be used to tailor 
    # the nucleon physics. If something other than the default CEM03 model is
    # desired for medium energy projecties, exactly one of the other models
    # should be switched on. In other words, only one of CEM, QMD, binary
    # cascade, INCL-ABLA cascade, or Bertini cascade should be set to True.
    # The others should be False. The standard Geant4 high precision (HP) 
    # neutrons can be switched off in favor the lower precision but faster
    # (LE) neutrons using the SetUseHPNeutrons(False) method. Note that this 
    # method replaces SetUseLENeutrons from earlier versions and reverses
    # the sense of the earlier switch. There is also a switch, which defaults
    # to True, to enforce a strict upper energy limit of 3GeV/u for the 
    # INCL-ABLA cascade. If this switch is False, the INCL-ABLA cascade model 
    # is used up to 8 GeV/u as are the Binary and Bertini cascades,
    # overlapping the QGSC model, which is applicable above 5 GeV/u. To 
    # fill the gap between the INCL-ABLA and QGSM models, the 
    # G4LEProtonInelastic model is used in the range 2.8-8 GeV/u, when 
    # the strict-energy-limit flag is set to True. Of course, all of these 
    # switches must be set to the desired position before mred.init().
    # After initialization, the values of these parameters are never used.
    if not options.disableNuclear:
        mred.physics.addModule('NucleonHadronic')
        if not options.particle == "neutron":
            mred.physics.module_dict['NucleonHadronic'].SetUseHPNeutrons(False)
    # mred.physics.module_dict['NucleonHadronic'].SetUseCEM(True)
    # mred.physics.module_dict['NucleonHadronic'].SetProtonsUseQMD(False)
    # mred.physics.module_dict['NucleonHadronic'].SetUseBinaryCascade(False)
    # mred.physics.module_dict['NucleonHadronic'].SetUseInclAblaCascade(False)
    # mred.physics.module_dict['NucleonHadronic'].SetInclAblaUseStrictEnergyLimit(True)
    # mred.physics.module_dict['NucleonHadronic'].SetUseBertiniCascade(False)
    # mred.physics.module_dict['NucleonHadronic'].SetUseHPNeutrons(True)
    #
    # mred.physics.addModule('NucleonHadronicJQMD')	# Experimental, Deprecated
    # mred.physics.addModule('NucleonHadronicPHITS')	# Experimental, Deprecated
    #
    # Deprecated "NucleonInelasticA" = QGSP + binary cascade + HP neutrons
    # Deprecated "NucleonInelasticB" = QGSP + Bertini cascade + HP neutrons
    # Deprecated "NucleonInelasticC" = QGSP + Bertini cascade + standard neutrons
    # Deprecated "NucleonInelasticD" = QGSP + binary cascade + standard neutrons
    #
    # mred.physics.addModule('NucleonInelasticA')	# Deprecated
    # mred.physics.addModule('NucleonInelasticB')	# Deprecated
    # mred.physics.addModule('NucleonInelasticC')	# Deprecated
    # mred.physics.addModule('NucleonInelasticD')	# Deprecated
    # mred.physics.addModule('NucleonInelasticIA') # Experimental INCL-ABLA module
    # --------------------------------------------------------------------------
    # 'IonInelasticLAQGSM' is the default ion-ion model
    #
    if not options.disableNuclear:
        mred.physics.addModule('IonInelasticLAQGSM')  # This is the recommended ion module.
    # mred.physics.addModule('IonInelasticDPMJET') # This is for the web site.
    # mred.physics.addModule('IonInelastic')       # This uses LAQGSM but also with QMD.
    #
    # The following are for using G4QMD
    # mred.physics.addModule('IonInelasticG4QMD')
    # The following can be used together to pick the FRAG model to replace GEM.
    # This is a one-way choice. It is not possible to go back to GEM.
    # mred.physics.module_dict['IonInelasticG4QMD'].GetG4QMDModel().UnUseGEM()
    # mred.physics.module_dict['IonInelasticG4QMD'].GetG4QMDModel().UseFRAG()
    #
    # These are experimental and not for general use.
    # mred.physics.addModule('IonPhysics')	# G4 module. Experimental. No inelastic scattering.
    # mred.physics.addModule('IonQMD')		# G4 module. Experimental. Max 10GeV/u.
    #
    # These are deprecated.
    # mred.physics.addModule('AltIonInelastic')
    # --------------------------------------------------------------------------
    # The 'Decay' module is automatically included when processes are constructed.
    mred.physics.module_dict['Decay'].SetIncludeRadioactiveDecay(False)
    # --------------------------------------------------------------------------
    # These can be used to turn on and set the threshold for Penelope transport.
    # mred.physics.module_dict['EmStandardScreened'].SetUseFortranPenelope(True)
    # mred.physics.module_dict['EmStandardScreened'].SetPenelopeThreshold(50.*keV)
    # Or to use Penelope2008 and Geant4 Penelope together...
    # mred.physics.module_dict['EmPenelopeQED'].SetUseFortranPenelope(True)
    # mred.physics.module_dict['EmPenelopeQED'].SetPenelopeThreshold(50.*keV)
    # --------------------------------------------------------------------------
  
# ------------------------------------------------------------------------------
mred.physics.use_multiplier_primary_only=True
mred.physics.sigma_multiplier=options.sigmaBias
mred.physics.use_track_weighting=True
print "cross_section_wrapper_info = %s" % mred.physics.cross_section_wrapper_info
#print mred.physics.list.GetSigmaBiasInfo()

# -----------------------------------------------------------------------------
mred.physics.range_cuts = options.rangeCuts*micrometer
# -----------------------------------------------------------------------------


###############################################################################
# Create a layered device
###############################################################################
d=mred.setDevice('rpp')
mred.materials.enableBasicElectronicMaterials()
#mred.materials.enableMaterial('air')
#rmed.materials.enableMaterial('lead')

# Create the layered device for an array of devices. This stack is created to
# include extra material around the array.  Because of this, the sensitive
# detector is not the entire width of the silicon layer. I have also omitted
# the oxide material for STI instead lumping it into the silicon layer. This
# will hopefully avoid the boundary stepping problem with energy deposition
# seen in other simulations.
layers = []
nCols=options.nCols
nRows=options.nRows

margin=0*micrometer
substrateDepth=5*micrometer

# The dimension of the cell and sensitive volume are based on shrink from
# a 65nm SRAM cell. Here we calculate the dimensions needed for the
# geometric construction. The gate width and length are assumed to be the same
# as the feature size.
scaleFactor = options.tech / 65.0
cellLength=475*nm * scaleFactor
cellWidth=975*nm * scaleFactor
gateLength=options.tech*nm
gateWidth=1*gateLength

if options.process == 'bulk':
    # In bulk, the charge collection volume is related to the gate
    # length, gate width, and assumed to be 0.5 um deep due to the well
    # structure.
    svLength=256*nm * scaleFactor # 3*gateLength # [S-G-D]
    svWidth=115*nm * scaleFactor # gateWidth
    svDepth=0.5*micrometer
elif options.process == 'soi':
    # In SOI, the charge collection volume is related to the gate
    # length, gate width, and assumed to be 10 nm deep due to the isolation
    svLength=256*nm * scaleFactor # 3*gateLength # [S-G-D]
    svWidth=115*nm * scaleFactor # gateWidth
    svDepth=0.010*micrometer
elif options.process == 'bulkFinFET':
    # In bulk, the charge collection volume is related to the gate
    # length, gate width, and assumed to be 0.5 um deep due to the well
    # structure.
    svLength=256*nm * scaleFactor # 3*gateLength # [S-G-D]
    svWidth=115*nm * scaleFactor # gateWidth
    svDepth=0.240*micrometer
else:
    raise("Unknown process")

###################################################################################
# Solid Model Construction
###################################################################################
layerWidth=(nCols * 0.5) * micrometer + margin
layerLength=(nRows * 0.5) * micrometer + margin
print "Cell Size %g x %g x %g" % ( svWidth, svLength, svDepth )

layers.append(((layerWidth, layerLength, 200*micrometer), 'silicon')) # Enough silicon to stop 1MeV muons
layers.append(((layerWidth, layerLength, svDepth), 'silicon', 'default'))
layers.append(((layerWidth, layerLength, substrateDepth), 'silicon'))
d.setLayers(layers)
d.wafer_material='silicon'
mred.device.register(d.g4PV())

mred.init()
mred.beamOn(0)
#mred.physics.setProduceAugerElectrons(True)
#mred.physics.setProduceFluorescencePhotons(True)
#if options.disableRecoils:
#    mred.physics.screened_scattering.allow_energy_deposition = False
# Enable event energy histogramming for the SD.  Here we create a custom
# histogram instead of the default associated with the detector. The native
# sensitive volume model is insufficient to model multiple devices, unless
# there are multiple detectors. To handle this, a common histogram is created
# and single event mode will be responsible for tallies based on the energy
# deposited in each device.
sd = mred.sd_vector[0]
print "Configuring sensitive detector %s" % ( sd.name )
mred.runAct.SetAutoGenerateHistogramsFlag(False)
mred.accumulate_histograms=False

commonHistogram=mred.getNewHistogram(hMin=0.1*keV,hMax=100*MeV,nBins=600,logSpacing=True)
commonHistogram.clear()

# Create each device as a group of associated sensitive volumes
print "Wafer dimensions: %s" % (str(d.wafer_dimensions))
zoffset=(d.wafer_dimensions[2]/ 2.0) - (substrateDepth+svDepth)

for row in range(nRows):
    for col in range(nCols):
        #device=Device('r%d_c%d' % (row, col))
        #device.group=sd.addSensitiveVolume('group')
        #device.group.on=True
        #device.group.weight=1
        #device.Ecrit=device.group.min_valid_total_energy=options.min_energy*MeV

        xoffset=(col - nCols/2.0)*cellLength
        yoffset=(row - nRows/2.0)*cellWidth
        
        sv=sd.addSensitiveVolume('rpp')
        sv.setCenter(vector3d(xoffset,yoffset,zoffset+svDepth/2.0))
        sv.setSize(vector3d(svLength,svWidth,svDepth))
        sv.on=True
        sv.weight=1
        sv.cpp_sv.SetDisplayColour(PyG4Core.G4Colour(0,1,0,min(1,sv.weight)))
        deviceList.append(sv)
        #print "Created SV size %s at %s" % (vector3d(svLength,svWidth,svDepth),vector3d(xoffset,yoffset,zoffset+svDepth/2.0), )

# This is really only enabled for efficiency. We don't want to enter the single
# event callback for every particle, so putting a small, but meaningful filter
# here will help a lot. If the model has been calibrated, this should be set
# just below the critical energy.
sd.use_weighted_as_total=True
sd.min_valid_total_energy=options.min_energy*MeV

mred.runMgr.SetVerboseLevel(1)

################################################################################
# HDF5 Output
################################################################################
if batch_vars:
    mred.hdf5.write_output_files=True
    mred.hdf5.include_histograms=True
    mred.hdf5.include_tracks=options.exportTracks
    mred.hdf5.file_path='hdf5_output/%s' % (batch_vars.runName)
    mred.hdf5.file_name="%s.%s.%03d.hdf5" % (batch_vars.runName, batch_vars.isotime, batch_vars.index)
    pickleFile="%s.%s.%03d.pickle" % (batch_vars.runName, batch_vars.isotime, batch_vars.index)
elif options.runName:
    mred.hdf5.write_output_files=True
    mred.hdf5.include_histograms=True
    mred.hdf5.include_tracks=options.exportTracks
    mred.hdf5.file_path='hdf5_output'
    mred.hdf5.file_name="%s.hdf5" % (options.runName)
    pickleFile="%s.pickle" % (options.runName)

################################################################################
# Set up the gun
################################################################################
if options.environment:
    mred.gun.random_spatial_sampling='isotropicFlux'
    mred.gun.random_energy_sampling='log'
elif options.spectrum:
    mred.gun.random_spatial_sampling='directionalFlux'
    mred.gun.energy_spectrum_file='%s' % (options.spectrum)
    mred.gun.random_energy_sampling='log'
    if options.particle:
        mred.gun.setParticle(options.particle)
else:
    if options.particle == 'ion' and options.beamZ == 1 and options.beamA == 1:
        mred.physics.screened_scattering.mean_free_path_scale_factor=0.1
        mred.gun.setParticle('proton')
    elif options.particle == 'ion':
        mred.gun.setParticle('ion', options.beamZ, options.beamA)
    else:
        mred.gun.setParticle(options.particle)
    if options.beamP:
        M=mred.gun.particle.GetPDGMass()
        options.beamE = -M + math.sqrt( M**2 + (options.beamP)**2 )
    mred.gun.energy=options.beamE*MeV
    mred.gun.random_spatial_sampling='directionalFlux'

if mred.gun.random_spatial_sampling=='directionalFlux':
    # Set the gun direction
    mred.gun.random_use_device_radius=False
    mred.gun.setRandomBox()

def singleEventCallback(evt):
    """This function will serve to evaluate the current event for each device."""
    # Get the energy deposited and event weight
    evtId=mred.runMgr.GetEventCount()-1
    evtWeight=mred.evtAct.ComputeEventWeight(evt)

    for device in deviceList:
        Edep=device.ionizing_energy
        commonHistogram.add(Edep,evtWeight/float(nRows*nCols))

    # Export track and random seeds for event
    if options.exportTracks:                
        idx=commonHistogram.index(Edep)
        counts=commonHistogram.y()[idx]
        if counts <= 5.0/(nRows*nCols):
            seeds=mred.last_random_seeds
            print "E= %g (%d, %d)" % ( Edep, seeds[0], seeds[1])
            sys.stdout.flush()
            if mred.hdf5.write_output_files:
                mred.hdf5.writeEvent(evt)

    # Draw the event if in interactive mode
    if options.enableDX:
        mred.dx.displayMredEvent(evt)
    if options.enableSUV:
        mred.evtAct.DrawEvent(evt)

################################################################################
# Viewers
################################################################################
if options.enableSUV:
    mred.suv()
if options.enableDX or options.exportDX:
    mred.dx.captureGeometry([mred.detCon.GetPhysicalWorld(),sd.g4PV()],custom_color_map=dx_color_dict,opacity_multiplier=1)

################################################################################
# Filters
################################################################################
if options.particle == "neutron":
    mred.filter.on=True
    mred.filter.verbose=False
    mred.filter.include_particles=True
    mred.filter.include_scatters=True
    mred.filter.include_nuclear_elastic=True
    mred.filter.include_coulomb_elastic=True

################################################################################
# Run the simulation
################################################################################
def run(nIons):
    """A very convenient function for running in interactive mode"""
    callbackFn=singleEventCallback
    # Prevent the process from being killed if it exceeds the requested run time
    # Instead, abort the simulation and exit gracefully.
    if batch_vars:
        mred.runMgr.SetMaxRunSeconds(int(batch_vars.runTime)-60)
    mred.progress_interval=int(nIons/10.0)

    if options.environment:
        mred.runRadiationEnvironment(nIons,options.environment,function=callbackFn)
    else:
        mred.runSingleEventMode(nIons,callbackFn)

def replay(seeds):
    """A very convenient function for replaying events"""
    mred.random_seeds=seeds
    mred.filter.verbose=True
    run(1)
    mred.filter.verbose=False

###############################################################################
# Main
###############################################################################
if options.random_seeds:
    print "Replaying event %s" % (options.random_seeds)
    #replay(eval(options.random_seeds))
    mred.random_seeds=eval(options.random_seeds)
    run(options.nIons)
else:
    run(options.nIons)

PostProcess()
print "Hello?"

#!/usr/bin/env python3
# (C) 2017 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.

#############################################################################
# Depicts the unpaired and clash interactions of an active site
#############################################################################

'''
Usage: python unpairedmap2img.py -complex production.pdb -out (lig).svg
'''

import sys
from openeye import oechem
from openeye import oedepict
from openeye import oegrapheme
from openeye.oespicoli import *
oechem.OESetLicenseFile("/home/victor/oe_license.txt")

def main(argv=[__name__]):

    itf = oechem.OEInterface()
    oechem.OEConfigure(itf, InterfaceData)
    oedepict.OEConfigureImageWidth(itf, 900.0)
    oedepict.OEConfigureImageHeight(itf, 600.0)
    oedepict.OEConfigure2DMolDisplayOptions(itf, oedepict.OE2DMolDisplaySetup_AromaticStyle)
    oechem.OEConfigureSplitMolComplexOptions(itf, oechem.OESplitMolComplexSetup_LigName |
                                             oechem.OESplitMolComplexSetup_CovLig)

    if not oechem.OEParseCommandLine(itf, argv):
        return 1

    if itf.HasString("-complex") and (itf.HasString("-protein") or itf.HasString("-ligand")):
        oechem.OEThrow.Warning("Only complex in %s file fill be used!" % itf.GetString("-complex"))

    if not (itf.HasString("-complex")) ^ (itf.HasString("-protein") and itf.HasString("-ligand")):
        oechem.OEThrow.Fatal("Please specify either complex or ligand and protein input files!")

    oname = itf.GetString("-out")

    ext = oechem.OEGetFileExtension(oname)
    if not oedepict.OEIsRegisteredImageFile(ext):
        oechem.OEThrow.Fatal("Unknown image type!")

    ofs = oechem.oeofstream()
    if not ofs.open(oname):
        oechem.OEThrow.Fatal("Cannot open output file!")

    # initialize protein and ligand

    protein = oechem.OEGraphMol()
    ligand = oechem.OEGraphMol()
    if not get_protein_and_ligand(protein, ligand, itf):
        oechem.OEThrow.Fatal("Cannot initialize protein and/or ligand!")

    # depict unpaired interaction map

    width, height = oedepict.OEGetImageWidth(itf), oedepict.OEGetImageHeight(itf)
    image = oedepict.OEImage(width, height)

    interactive_legend = itf.GetBool("-interactive-legend")

    cwidth, cheight = width, height
    if not interactive_legend:
        cwidth = cwidth * 0.8

    opts = oegrapheme.OE2DActiveSiteDisplayOptions(cwidth, cheight)
    opts.SetRenderInteractiveLegend(interactive_legend)
    oedepict.OESetup2DMolDisplayOptions(opts, itf)

    if interactive_legend:
        depict_unpairedmap(image, protein, ligand, opts)
    else:
        main_frame = oedepict.OEImageFrame(image, width * 0.80, height,
                                           oedepict.OE2DPoint(width * 0.2, 0.0))
        legend_frame = oedepict.OEImageFrame(image, width * 0.20, height,
                                             oedepict.OE2DPoint(width * 0.0, 0.0))
        depict_unpairedmap(main_frame, protein, ligand, opts, legend_frame)

    if ext == 'svg' and interactive_legend:
        iconscale = 0.5
        oedepict.OEAddInteractiveIcon(image, oedepict.OEIconLocation_TopRight, iconscale)
    oedepict.OEDrawCurvedBorder(image, oedepict.OELightGreyPen, 10.0)

    oedepict.OEWriteImage(oname, image)

    return 0


def depict_unpairedmap(image, protein, ligand, opts, legend_frame=None):
    """
    :type image: oedepict.OEImageBase
    :type protein: oechem.OEMolBase
    :type ligand: oechem.OEMolBase
    :type opts: oedepict.OE2DMolDisplayOptions
    :type legend_frame: oedepict.OEImageBase
    """

    # perceive interactions

    asite = oechem.OEInteractionHintContainer(protein, ligand)
    if not asite.IsValid():
        oechem.OEThrow.Fatal("Cannot initialize active site!")
    asite.SetTitle(ligand.GetTitle())

    oechem.OEPerceiveInteractionHints(asite)

    # depiction

    oegrapheme.OEPrepareActiveSiteDepiction(asite)
    adisp = oegrapheme.OE2DActiveSiteDisplay(asite, opts)
    oegrapheme.OERenderUnpairedInteractionMap(image, adisp)

    if legend_frame is not None:
        lopts = oegrapheme.OE2DActiveSiteLegendDisplayOptions(12, 1)
        oegrapheme.OEDrawUnpairedInteractionMapLegend(legend_frame, adisp, lopts)


def split_complex(protein, ligand, sopts, complexmol):

    water = oechem.OEGraphMol()
    other = oechem.OEGraphMol()

    pfilter = sopts.GetProteinFilter()
    wfilter = sopts.GetWaterFilter()
    sopts.SetProteinFilter(oechem.OEOrRoleSet(pfilter, wfilter))
    filtercategory = oechem.OEMolComplexFilterCategory_Nothing
    sopts.SetWaterFilter(oechem.OEMolComplexFilterFactory(filtercategory))

    oechem.OESplitMolComplex(ligand, protein, water, other, complexmol, sopts)

    return ligand.NumAtoms() != 0 and protein.NumAtoms() != 0


def get_protein_and_ligand(protein, ligand, itf):

    if itf.HasString("-complex"):

        # separate ligand and protein in complex

        iname = itf.GetString("-complex")

        ifs = oechem.oemolistream()
        if not ifs.open(iname):
            oechem.OEThrow.Fatal("Cannot open input complex file!")

        complexmol = oechem.OEGraphMol()
        if not oechem.OEReadMolecule(ifs, complexmol):
            oechem.OEThrow.Fatal("Unable to read complex from %s" % iname)

        if not oechem.OEHasResidues(complexmol):
            oechem.OEPerceiveResidues(complexmol, oechem.OEPreserveResInfo_All)

        sopts = oechem.OESplitMolComplexOptions()
        oechem.OESetupSplitMolComplexOptions(sopts, itf)

        if not split_complex(protein, ligand, sopts, complexmol):
            oechem.OEThrow.Fatal("Cannot separate complex!")
    else:

        # read ligand and protein from separate files

        pname = itf.GetString("-protein")

        ifs = oechem.oemolistream()
        if not ifs.open(pname):
            oechem.OEThrow.Fatal("Cannot open input protein file!")

        if not oechem.OEReadMolecule(ifs, protein):
            oechem.OEThrow.Fatal("Unable to read protein from %s" % pname)

        lname = itf.GetString("-ligand")

        ifs = oechem.oemolistream()
        if not ifs.open(lname):
            oechem.OEThrow.Fatal("Cannot open input ligand file!")

        if not oechem.OEReadMolecule(ifs, ligand):
            oechem.OEThrow.Fatal("Unable to read ligand from %s" % lname)

    return ligand.NumAtoms() != 0 and protein.NumAtoms() != 0


#############################################################################
# INTERFACE
#############################################################################
InterfaceData = '''

!CATEGORY "input/output options :" 1

  !PARAMETER -complex 1
    !ALIAS -c
    !TYPE string
    !REQUIRED false
    !VISIBILITY simple
    !BRIEF Input filename of the protein-ligand complex
  !END

  !PARAMETER -protein 2
    !ALIAS -p
    !TYPE string
    !REQUIRED false
    !VISIBILITY simple
    !BRIEF Input filename of the protein
  !END

  !PARAMETER -ligand 3
    !ALIAS -l
    !TYPE string
    !REQUIRED false
    !VISIBILITY simple
    !BRIEF Input filename of the ligand
  !END

  !PARAMETER -out 4
    !ALIAS -o
    !TYPE string
    !REQUIRED true
    !VISIBILITY simple
    !BRIEF Output filename of the generated image
  !END

!END

!CATEGORY "active site display options:" 2

  !PARAMETER -interactive-legend
    !ALIAS -ilegend
    !TYPE bool
    !DEFAULT false
    !REQUIRED false
    !VISIBILITY simple
    !BRIEF Visualize legend on mouse hover (SVG feature)
  !END

!END
'''

if __name__ == "__main__":
    sys.exit(main(sys.argv))

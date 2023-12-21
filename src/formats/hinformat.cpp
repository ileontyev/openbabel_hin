/**********************************************************************
Copyright (C) 2000 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>

#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>

#include <cstdlib>

#include <iostream>// IL Modified
#include <openbabel/generic.h>  // IL Modified

using namespace std;
  char* strcat_mem(const char *str1, const char *str2) {
   char* result;
   //cout << "START strcat_mem"  << endl;
    int l2 = strlen(str2) ;
    if (!str1) {
     result = (char*) malloc(l2+1);
     memcpy(result, str2, l2) ;
     result[l2] = '\0' ;
     return result;
    }

    int l1 = strlen(str1) ;

    result = (char*) malloc(l1 + l2 + 1);
    if(!result) return result;
    //cout << "memcpy str1:" << str1 << endl;
    memcpy(result, str1, l1) ;
    //cout << "memcpy str2:" << str2 << endl;
    memcpy(result + l1, str2, l2 );
    result[l1 + l2] = '\0' ;
    //cout << "return strcat_mem:" << result << endl;
    //*str1  = *result;
    return result;
  }

namespace OpenBabel
{

  class HINFormat : public OBMoleculeFormat
  {
  public:
   // IL modifications
    //Register this format type ID
    HINFormat()
    {
      OBConversion::RegisterFormat("hin",this, "chemical/x-hin");
      OBConversion::RegisterOptionParam("r", this, 0, OBConversion::INOPTIONS);
      OBConversion::RegisterOptionParam("T", this, 0, OBConversion::INOPTIONS); // "t"-option is already registered in fastsearchformat.cpp
      OBConversion::RegisterOptionParam("c", this, 0, OBConversion::INOPTIONS);
      OBConversion::RegisterOptionParam("c", this, 0, OBConversion::OUTOPTIONS);
      OBConversion::RegisterOptionParam("l", this, 0, OBConversion::OUTOPTIONS);
      OBConversion::RegisterOptionParam("u", this, 0, OBConversion::OUTOPTIONS);
    }
    const char* Description() override  // required
    {
      return "HyperChem HIN format\n"
        "Read Options e.g. -aR -aT\n"
        "  R               Persive residue information (atom label, resname, chainID) automatically, ignoring info parsed from HIN. Works for PROTEIN only.\n"
        "  T               Read atom Types from HIN (glitchy, because convertion formats usually Translate() types using predifined ttab).\n"
        "                  By default, atom types from HIN are ignored, then auto persived at the converted format Write routine (e.g. to mol2 types)\n"
        "  C               Parse pair-formatted comments (fmt: ;field-name value) to OBPairData type, otherwise multi-line Generic OBCommentData\n\n"
        "Write Options e.g. -xl\n"
        "  l               (Placeholder) Output ignores residue information (only ligands)\n"
        "  c               (Placeholder) Write UCSF Dock scores saved in comments preceding molecules\n"
        "  u               (Placeholder) Do not write formal charge information in UNITY records\n\n";
    }
    const char* SpecificationURL() override
    { return "\nHyperChem7 Manual - Appendix-D HIN files: http://www.chemistry-software.com/pdf/Hyperchem_full_manual.pdf\nOr at\nhttps://wiki.jmol.org/index.php/File_formats/Formats/HIN"; }  // optional
 // End of IL modifications





    const char* GetMIMEType() override
    { return "chemical/x-hin"; }  // optional

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    bool ReadMolecule(OBBase* pOb, OBConversion* pConv) override;
    bool WriteMolecule(OBBase* pOb, OBConversion* pConv) override;
  };
  //***
  
  //Make an instance of the format class
  HINFormat theHINFormat;

  /////////////////////////////////////////////////////////////////
  bool HINFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if (pmol == nullptr)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    // Right now only read in the first molecule
    int i;
    int max, bo;
    char buffer[BUFF_SIZE];
    string str,str1;
    double x,y,z;
    OBAtom *atom;
    vector<string> vs,line;
    char *comment = nullptr; // IL Modified 
    int len;

    ifs.getline(buffer, BUFF_SIZE);
    while (ifs.good() && (strstr(buffer, "mol") == nullptr || buffer[0] == ';')) //The "mol" in comment line should be ignored.
      {
    // IL modifications - Push comments to mol.SetData using OBCommentData *cd
        if ( buffer[0] == ';' )
          {
            if (comment) {
              comment = strcat_mem(comment, "\n");
            }
            comment = strcat_mem(comment, buffer + 1); // buffer + 1 - skip buffer[0](";") : https://stackoverflow.com/questions/5711490/c-remove-the-first-character-of-an-array
          }    // End of IL modifications
    
        ifs.getline(buffer, BUFF_SIZE);
        if (ifs.peek() == EOF || !ifs.good())
          return false;
      }

    // IL modifications
    tokenize(vs,buffer);

    //printf("Address of pointer *mol is %p\n", mol);
    //printf("buffer = %s, vs[2] = \"%s\"\n", buffer,(char*)vs[2].c_str());

    char *mol_title = nullptr;
    if ( strcmp((char*)vs[0].c_str(), "mol") == 0 ) {
      if (vs.size() == 3){    // Don't really know how long it'll be
        mol_title = strcat_mem(mol_title, (char*)vs[2].c_str());
        //cout << "vs[2]: " << vs[2] <<" mol_title: " << mol_title << " title: " << title << endl;
        //printf("Mol line: mol_title = %s, title = %s\n", mol_title, title);
      }
    }    // End of IL modifications

    ifs.getline(buffer, BUFF_SIZE);
    if (!ifs.good())
      return false; // ended early

    // We need to prevent chains perception routines from running while
    // we are adding residues from the PDB file
    mol.SetChainsPerceived();

    mol.BeginModify();

    while (ifs.good() && strstr(buffer, "endmol") == nullptr)
      {
        if(buffer[0]==';'){
          // IL modifications - Push comments to mol.SetData using OBCommentData *cd
          if (comment) {
            comment = strcat_mem(comment, "\n");
          }
          comment = strcat_mem(comment, buffer + 1); // buffer + 1 - skip buffer[0](";") : https://stackoverflow.com/questions/5711490/c-remove-the-first-character-of-an-array
          // End of IL modifications


           ifs.getline(buffer, BUFF_SIZE);
           continue; //The comment Line in HIN should be ignored.
        }

        tokenize(vs,buffer);

        // IL modifications
        // HIN Format Specification is in HyperChem7 manual, AppendixD: https://wiki.jmol.org/index.php/File_formats/Formats/HIN
        if (vs.size() < 3)    // Don't really know how long it'll be
          {
            ifs.getline(buffer, BUFF_SIZE);
            continue;
          }
        OBResidue *res  = (mol.NumResidues() > 0) ? mol.GetResidue(mol.NumResidues()-1) : nullptr;
        if ( strcmp((char*)vs[0].c_str(), "res") == 0 ) {
          string resnum = (char *)vs[1].c_str();
          string resname = (char *)vs[2].c_str();
          res = mol.NewResidue();
          res->SetName(resname);
          res->SetNum(resnum);
          // res->SetInsertionCode(insertioncode);
          char chain;
          if (vs.size() >= 6) {
            chain = (char)vs[5].c_str()[0];
            res->SetChain(chain);
          }
          /*
          cout << buffer << endl;
          cout << "Residue-" << resname << resnum << ", chain-" << chain << endl;
          cout << "SetResidue-" << res->GetNum() << ": " << (char *)res->GetName().c_str() << ", chain-" << res->GetChain() << endl;
          printf("SetResidue-%d: %s, chain-%c\n",res->GetNum(),(char*)res->GetName().c_str(),(char)res->GetChain());
          printf("Residue-%s: %s, chain-%c\n",resnum.c_str(),resname.c_str(), chain);
          */

          ifs.getline(buffer, BUFF_SIZE);
          continue;

        }
        // End of IL modifications

        if (vs.size() < 11)    // Don't really know how long it'll be
          {
            ifs.getline(buffer, BUFF_SIZE);
            continue;
          }

        atom = mol.NewAtom();
        atom->SetAtomicNum(OBElements::GetAtomicNum(vs[3].c_str()));
        atom->SetPartialCharge(atof(vs[6].c_str()));
        x = atof((char*)vs[7].c_str());
        y = atof((char*)vs[8].c_str());
        z = atof((char*)vs[9].c_str());
        atom->SetVector(x,y,z);

        max = 11 + 2 * atoi((char *)vs[10].c_str());
        for (i = 11; i < max; i+=2)
          {
            switch(((char*)vs[i+1].c_str())[0]) // First char in next token
              {
              case 's':
                bo = 1;
                break;
              case 'd':
                bo = 2;
                break;
              case 't':
                bo = 3;
                break;
              case 'a':
                bo = 5;
                break;
              default :
                bo = 1;
                break;
              }
            mol.AddBond(mol.NumAtoms(), atoi((char *)vs[i].c_str()), bo);
          }

     // IL modifications
        atom->SetFormalCharge(atof(vs[6].c_str())); // Unlike HyperChem, InterX uses vs[6] for a formal charge not a patial charge as set above (atom->SetPartialCharge)
        if (pConv->IsOption("T", OBConversion::INOPTIONS) != nullptr){
          atom->SetType((char *)vs[4].c_str());;   // Set Types only if auto preception is not forced by "T" option, otherwise ttab.Translate is used, like in MOL2.
          // Description of the atom Type translation ttab.Translate: https://openbabel.github.io/api/3.0/classOpenBabel_1_1OBTypeTable.shtml#a9b40b242258cc488f341d06412a1fb89
        }

        // Add Default Residue Info if it was not defined in HIN
        //OBAtom *atom = mol.GetAtom(mol.NumAtoms());
        if (res == nullptr) {
          res = mol.NewResidue();
          res->SetName("UNL");
          res->SetNum("1");
          // res->SetInsertionCode(insertioncode);
        }
        res->AddAtom(atom);
        res->SetSerialNum(atom, atoi((char*)vs[1].c_str()));
        res->SetAtomID(atom, (char *)vs[2].c_str());
/*
        OBAtom *atom1 = mol.GetAtom(i);
        printf("Atom-%d: %s %s\n",res->GetSerialNum(atom), (char*)res->GetAtomID(atom).c_str(), atom->GetType());
        printf("Residue-%s: %s\n",res->GetNumString(),(char*)res->GetName().c_str());
*/
    // End of IL modifications
          
        ifs.getline(buffer, BUFF_SIZE);
      }

    // clean out remaining blank lines
    // blank lines cleaning codes rewritten for avoiding peek() and tellg() bugs
    // https://github.com/openbabel/openbabel/issues/1569
    std::streampos ipos;
    do
    {
      ipos = ifs.tellg();
      ifs.getline(buffer,BUFF_SIZE);
    }
    while(strlen(buffer) == 0 && !ifs.eof() );
    ifs.seekg(ipos);

    mol.EndModify();

    mol.SetTitle(title); // Use filename as Title
    if (mol_title) {
      mol.SetTitle(mol_title); // Use molecule name parsed from HIN
    }

    mol.SetPartialChargesPerceived();
    //if (has_residue_information) - not used here, because SetChainsPerceived() is the only way to preserve AtomID (labels) info
    if (pConv->IsOption("R", OBConversion::INOPTIONS) == nullptr){
        // Mark as percieved (will use residue info from HIN), unless forced by the "r" option,
        // which ignores res info (res->SetAtomID(atom), resname, chainID) from input and try to auto persive. Auto percieve works for protein only.
      mol.SetChainsPerceived();
    }
    if (pConv->IsOption("T", OBConversion::INOPTIONS) != nullptr){
      mol.SetAtomTypesPerceived();   // Mark as percieved (will use types parsed from HIN), unless forced by the "t" option which will auto pecieve types.
    }

    // Parse comments to OBPairData type if read option -C or to Generic OBCommentData
    if (comment)
      {
        if (pConv->IsOption("C", OBConversion::INOPTIONS) != nullptr)
          {
            tokenize(line, comment,"\n");
            comment[0] = '\0'; // Empty the char array by defining its 1st element as end of string
            for (i = 0; i < line.size(); i+=1)
            {
              tokenize(vs,line[i]);
                  if (vs.size() == 2)
                  {
                     //char attr[32], val[32];
                     //sscanf(buffer, "########## %[^:]:%s", attr, val);
                     OBPairData *pd = new OBPairData;
                     pd->SetAttribute(vs[0]);
                     pd->SetValue(vs[1]);
                     pd->SetOrigin(fileformatInput);
                     //pd->SetOrigin(local);
                     mol.SetData(pd);
                   }else{
                     if (comment) {
                       comment = strcat_mem(comment, "\n");
                      }
                      comment = strcat_mem(comment, (char *)line[i].c_str() ); // line[i].c_str() - convert string to (const char*) type
                   }
            }
          } 
        // Use the generic single multi-line comment 
        //must add generic data after end modify - otherwise it will be blown away
        OBCommentData *cd = new OBCommentData;
        cd->SetData(comment);
        cd->SetOrigin(fileformatInput);
        mol.SetData(cd);
       }   // End of IL modifications

    // Print Debugging Info
  /*
    cout << "mol_title: " << mol_title << " title: " << title << endl;
    cout << "SetTitle: " << mol.GetTitle() << endl;
    printf("\nAddress of pointer *mol is %p\n", mol);
    for (i = 1; i <= mol.NumAtoms(); i++)
      {
         atom = mol.GetAtom(i);
         OBResidue *res1 = atom->GetResidue();
         printf("Molecule: Title = %s, mol_title = %s\n", mol.GetTitle(), mol_title);
         printf("%d Atom-%d: %s X = %8.5f\n", i, res1->GetSerialNum(atom), (char*)res1->GetAtomID(atom).c_str()), atom->GetX();
         printf("Residue-%d: %s\n",res1->GetNum(),(char*)res1->GetName().c_str());
      }
      OBResidue *res2;
      FOR_ATOMS_OF_MOL(atom, mol)
      {
         res2 = atom->GetResidue();
         cout << atom->GetIdx() << " ";
         cout << atom->GetAtomicNum() << " ";
         cout << res2->GetSerialNum(atom) << " ";
         cout << (char*)res2->GetAtomID(atom).c_str() << " ";
         cout << atom->GetVector() << endl;
      }
      cout << "res2->GetAtoms(): " << endl;
      //cout << res2->GetAtoms() << endl;
*/
    return(true);
  }

  ////////////////////////////////////////////////////////////////

  bool HINFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if (pmol == nullptr)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    unsigned int i, file_num = 1;
    string str,str1;
    char buffer[BUFF_SIZE];
    OBAtom *atom;
    OBBond *bond;
    vector<OBBond*>::iterator j;
    char bond_char;

    // Dump CommentData
    if (mol.HasData(OBGenericDataType::CommentData))
    {
      OBCommentData *cd = (OBCommentData*)mol.GetData(OBGenericDataType::CommentData);
      ofs << cd->GetData() << endl;
    }
    
    // make sure to escape titles in double quotes
    // PR#1501694
    ofs << "mol " << file_num << " \"" << mol.GetTitle() << "\"\n";

    for(i = 1;i <= mol.NumAtoms(); i++)
      {
        atom = mol.GetAtom(i);
        OBResidue *res = atom->GetResidue();
        char* atype = atom->GetType();
        if (!atype){
          atype = (char*)"**";
        }
        snprintf(buffer, BUFF_SIZE, "atom %d %s %-3s %s  - %8.5f %8.5f  %8.5f  %8.5f %d ",
                i,
                (char*)res->GetAtomID(atom).c_str(),
                OBElements::GetSymbol(atom->GetAtomicNum()),
                atype,
                atom->GetPartialCharge(),
                atom->GetX(),
                atom->GetY(),
                atom->GetZ(),
                atom->GetExplicitDegree());
/*
        snprintf(buffer, BUFF_SIZE, "atom %d - %-3s **  - %8.5f %8.5f  %8.5f  %8.5f %d ",
                i,
                OBElements::GetSymbol(atom->GetAtomicNum()),
                atom->GetPartialCharge(),
                atom->GetX(),
                atom->GetY(),
                atom->GetZ(),
                atom->GetExplicitDegree());
*/
        ofs << buffer;
        for (bond = atom->BeginBond(j); bond; bond = atom->NextBond(j))
          {
            switch(bond->GetBondOrder())
              {
              case 1 :
                bond_char = 's';
                break;
              case 2 :
                bond_char = 'd';
                break;
              case 3 :
                bond_char = 't';
                break;
              case 5 :
                bond_char = 'a';
                break;
              default:
                bond_char = 's';
                break;
              }
            if (bond->IsAromatic())
              bond_char = 'a';

            snprintf(buffer,BUFF_SIZE, "%d %c ", (bond->GetNbrAtom(atom))->GetIdx(), bond_char);
            ofs << buffer;
          }
        ofs << endl;
      }
    ofs << "endmol " << file_num << endl;
    return(true);
  }

} //namespace OpenBabel

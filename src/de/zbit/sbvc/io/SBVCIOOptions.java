/*
 * $Id$
 * $URL$
 * ---------------------------------------------------------------------
 * This file is part of SBVC, the systems biology visualizer and
 * converter. This tools is able to read a plethora of systems biology
 * file formats and convert them to an internal data structure.
 * These files can then be visualized, either using a simple graph
 * (KEGG-style) or using the SBGN-PD layout and rendering constraints.
 * Some currently supported IO formats are SBML (+qual, +layout), KGML,
 * BioPAX, SBGN, etc. Please visit the project homepage at
 * <http://www.cogsys.cs.uni-tuebingen.de/software/SBVC> to obtain the
 * latest version of SBVC.
 *
 * Copyright (C) 2012-2013 by the University of Tuebingen, Germany.
 *
 * SBVC is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation. A copy of the license
 * agreement is provided in the file named "LICENSE.txt" included with
 * this software distribution and also available online as
 * <http://www.gnu.org/licenses/lgpl-3.0-standalone.html>.
 * ---------------------------------------------------------------------
 */
package de.zbit.sbvc.io;

import java.io.File;

import de.zbit.io.filefilter.SBFileFilter;
import de.zbit.kegg.io.KEGGtranslatorIOOptions.Format;
import de.zbit.util.prefs.KeyProvider;
import de.zbit.util.prefs.Option;
import de.zbit.util.prefs.OptionGroup;
import de.zbit.util.prefs.Range;

/**
 * Commandline options
 * 
 * 
 * @author Finja B&uuml;chel
 * @version $Rev$
 */
public interface SBVCIOOptions  extends KeyProvider{


  /*
   * Most important options: input, output and file format.
   */
  
  /**
   * Path and name of the source, KGML formatted, XML-file.
   */
  public static final Option<File> INPUT = new Option<File>("INPUT",
      File.class,
      "Path and name of the source, OWL formatted, XML-file.",
//      TODO: Create range that accepts multiple file filter
//      KGML, SBML, PID XML, BioPAX      
      new Range<File>(File.class, SBFileFilter.createAllFileFilter()), (short) 2, "-i" );
      //new File(System.getProperty("user.dir")));

  /**
   * Path and name, where the translated file should be put.
   */
  public static final Option<File> OUTPUT = new Option<File>("OUTPUT",
      File.class,
      "Path and name, where the translated file should be put.",
      (short) 2, "-o" );//, new File(System.getProperty("user.dir")));

  /**
   * Target file format for the translation.
   */
  public static final Option<Format> FORMAT = new Option<Format>("FORMAT",
      Format.class, "Target file format for the translation.",
      new Range<Format>(Format.class, Range.toRangeString(Format.class)),
      (short) 2, "-f", Format.SBML);
  
  public static final Option<Integer> SPECIES = new Option<Integer>("SPECIES",
      Integer.class, 
      "The target species for the analysis is defined by the taxonomy id. " +
      "For instance, homo sapiens has the taxonomy id 9606",
//      Option.buildRange(String.class,
      // this is not possible, because the list length of 978 is too long!
      // Species.getListOfNames(ApplicationMain.getListOfSpecies(),Species.SCIENTIFIC_NAME)),
//               "{\"hsa\", \"mmu\"}"),
      (short) 2, "-s ", null);
  
//  /**
//   * if the file should be splitted
//   */
//  public static final Option<Boolean> SPLIT_MODE = new Option<Boolean>("SPLIT_MODE",
//      Boolean.class, "If this option is set true and the input file consists of several " +
//      		"pathways, for each pathway a result file is created. These files are returned " +
//      		"in one ZIP file.",
//      (short) 2, "-s", false);

  /**
   * Define the default input/ output files and the default output format.
   */
  @SuppressWarnings("unchecked")
  public static final OptionGroup<Object> BASE_OPTIONS = new OptionGroup<Object>(
      "Base options",
      "Define the default input/ output files and the conversion option.",
      INPUT, OUTPUT, FORMAT, SPECIES);
  
}

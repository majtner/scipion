/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Roberto Marabini (roberto@mipg.upenn.edu)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
/*****************************************************************************/
/* Job Description Language for the grid                                     */
/*****************************************************************************/

#ifndef _XMIPP_JDL_HH
#define _XMIPP_JDL_HH

#include <data/xmipp_funcs.h>

#include <vector>

/**@defgroup JDLFiles JDL Files
   @ingroup InterfaceLibrary */
//@{
/** JDL Files.
    This is a class to write Job Description Language files for the grid.
    */
class JDLFile
{
public:
    /// Type. By default "Job"
    std::string jdl_type;

    /// Jobtype. By default "normal"
    std::string jdl_jobtype;

    /// Virtual organization. By default "biomed"
    std::string jdl_virtualorganization;

    /** Job rootname.
        stderr and stdout are dumped into job_root.stderr.txt and
        job_root.stdout.txt. The invoked program is job_root.sh*/
    std::string job_root;

    /** Running Script.
        The script job_root.sh is generated from the information of this
        set of lines. */
    std::vector <std::string> job_sh;

    /** Running programs.
        These are the programs that must be sent to be executed by the
        script. */
    std::vector <std::string> running_programs;

    /** Input files. */
    std::vector <std::string> input_files;

    /** Input selfiles. */
    std::vector <std::string> input_selfiles;

    /** Input wildfiles. */
    std::vector <std::string> input_wildfiles;

    /** Output files. */
    std::vector <std::string> output_files;

    /** Output selfiles. */
    std::vector <std::string> output_selfiles;

    /** Output wildfiles. */
    std::vector <std::string> output_wildfiles;

    /** Where to publish the input data.
        If not given, the Storage Element (SE) of the Computing
        Element (CE) with more free CPUs is selected. */
    std::string publisher_data_in;

    /** Where to publish the output data.
        If there is no data_out publisher, but there is a data_in publisher,
        then the data_in publisher is also the data out. If none is given,
        then baudelaire.cnb.uam.es is the data_out publisher, and the
        data_in publisher is taken as the Storage Element (SE) of the Computing
        Element (CE) with more free CPUs. */
    std::string publisher_data_out;

    /** Local output directory */
    std::string local_output_dir;
public:
    JDLFile()
    {
        clear();
    }
public:
    /** Write the content of this object to a file.
        The script job_root.sh and job_root.jdl are created.*/
    void write();

    /** Empties actual JDL structure. */
    void clear();

    /** Add file to the input. */
    void add_file_to_input_files(const std::string &str)
    {
        input_files.push_back(str);
    }

    /** Add file to the input selfiles.
        All files in this selfile will be packed and sent to the remote host.*/
    void add_file_to_input_selfiles(const std::string &str)
    {
        input_selfiles.push_back(str);
    }

    /** Add file to the input wildfiles.
        All files with this pattern will be packed and sent to the remote host.*/
    void add_file_to_input_wildfiles(const std::string &str)
    {
        input_wildfiles.push_back(str);
    }

    /** Add file to the output. */
    void add_file_to_output_files(const std::string &str)
    {
        output_files.push_back(str);
    }

    /** Add file to the output selfiles.
        All files in this selfile will be packed and sent from the remote host.*/
    void add_file_to_output_selfiles(const std::string &str)
    {
        output_selfiles.push_back(str);
    }

    /** Add file to the output wildfiles.
        All files with this pattern will be packed and sent from the remote host.*/
    void add_file_to_output_wildfiles(const std::string &str)
    {
        output_wildfiles.push_back(str);
    }

    /** Add command to the running script. */
    void add_command_to_script(const std::string &str)
    {
        job_sh.push_back(str);
    }

    /** Add program to the running programs. */
    void add_program_to_running_programs(const std::string &str)
    {
        running_programs.push_back(str);
    }
};
//@}
#endif

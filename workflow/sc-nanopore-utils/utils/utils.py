""" This module contains helper functions mostly taken from the pbio project.
They are copied here to avoid installing the entire pbio project.
"""

import os
import sys

import logging
logger = logging.getLogger(__name__)

# logging_utils

def add_logging_options(parser, default_log_file=""):
    """ This function add options for logging to an argument parser. In 
        particular, it adds options for logging to a file, stdout and stderr.
        In addition, it adds options for controlling the logging level of each
        of the loggers, and a general option for controlling all of the loggers.

        Args:
            parser (argparse.ArgumentParser): an argument parser

        Returns:
            None, but the parser has the additional options added
    """

    logging_options = parser.add_argument_group("logging options")

    default_log_file = ""
    logging_level_choices = ['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
    default_logging_level = 'WARNING'
    default_specific_logging_level = 'NOTSET'

    logging_options.add_argument('--log-file', help="This option specifies a file to "
        "which logging statements will be written (in addition to stdout and "
        "stderr, if specified)", default=default_log_file)
    logging_options.add_argument('--log-stdout', help="If this flag is present, then "
        "logging statements will be written to stdout (in addition to a file "
        "and stderr, if specified)", action='store_true')
    logging_options.add_argument('--no-log-stderr', help="Unless this flag is present, then "
        "logging statements will be written to stderr (in addition to a file "
        "and stdout, if specified)", action='store_true')

    logging_options.add_argument('--logging-level', help="If this value is specified, "
        "then it will be used for all logs", choices=logging_level_choices,
        default=default_logging_level)
    logging_options.add_argument('--file-logging-level', help="The logging level to be "
        "used for the log file, if specified. This option overrides "
        "--logging-level.", choices=logging_level_choices, 
        default=default_specific_logging_level)
    logging_options.add_argument('--stdout-logging-level', help="The logging level to be "
        "used for the stdout log, if specified. This option overrides "
        "--logging-level.", choices=logging_level_choices, 
        default=default_specific_logging_level)
    logging_options.add_argument('--stderr-logging-level', help="The logging level to be "
        "used for the stderr log, if specified. This option overrides "
        "--logging-level.", choices=logging_level_choices, 
        default=default_specific_logging_level)


def update_logging(args, logger=None, 
        format_str='%(levelname)-8s %(name)-8s %(asctime)s : %(message)s'):

    """ This function interprets the logging options in args. Presumably, these
        were added to an argument parser using add_logging_options.

    Parameters
    ----------
    args: argparse.Namespace
        a namespace with the arguments added by add_logging_options

    logger: logging.Logger or None
        a logger which will be updated. If None is given, then the default
        logger will be updated.

    format_str: string
        The logging format string. Please see the python logging documentation
        for examples and more description.

    Returns
    -------
    None, but the default (or given) logger is updated to take into account
        the specified logging options
    """

    # find the root logger if another logger is not specified
    if logger is None:
        logger = logging.getLogger('')
            
    logger.handlers = []

    # set the base logging level
    level = logging.getLevelName(args.logging_level)
    logger.setLevel(level)

    # now, check the specific loggers

    if len(args.log_file) > 0:
        h = logging.FileHandler(args.log_file)
        formatter = logging.Formatter(format_str)
        h.setFormatter(formatter)
        if args.file_logging_level != 'NOTSET':
            l = logging.getLevelName(args.file_logging_level)
            h.setLevel(l)
        logger.addHandler(h)

    if args.log_stdout:
        h = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter(format_str)
        h.setFormatter(formatter)
        if args.stdout_logging_level != 'NOTSET':
            l = logging.getLevelName(args.stdout_logging_level)
            h.setLevel(l)
        logger.addHandler(h)

    log_stderr = not args.no_log_stderr
    if log_stderr:
        h = logging.StreamHandler(sys.stderr)
        formatter = logging.Formatter(format_str)
        h.setFormatter(formatter)
        if args.stderr_logging_level != 'NOTSET':
            l = logging.getLevelName(args.stderr_logging_level)
            h.setLevel(l)
        logger.addHandler(h)
        
        
def get_logging_options_string(args):
    """ This function extracts the flags and options specified for logging options
        added with add_logging_options. Presumably, this is used in "process-all"
        scripts where we need to pass the logging options to the "process" script.

        Args:
            args (namespace): a namespace with the arguments added by add_logging_options

        Returns:
            string: a string containing all logging flags and options

    """

    args_dict = vars(args)

    # first, pull out the text arguments
    logging_options = ['log_file', 'logging_level', 'file_logging_level',
        'stdout_logging_level', 'stderr_logging_level']

    # create a new dictionary mapping from the flag to the value
    logging_flags_and_vals = {'--{}'.format(o.replace('_', '-')) : args_dict[o] 
        for o in logging_options if len(args_dict[o]) > 0}

    s = ' '.join("{} {}".format(k,v) for k,v in logging_flags_and_vals.items())

    # and check the flags
    if args.log_stdout:
        s = "--log-stdout {}".format(s)

    if args.no_log_stderr:
        s = "--no-log-stderr {}".format(s)

    return s


# shell_utils

import subprocess


def check_programs_exist(programs, raise_on_error=True, package_name=None, 
            logger=logger):

    """ This function checks that all of the programs in the list cam be
        called from python. After checking all of the programs, an exception
        is raised if any of them are not callable. Optionally, only a warning
        is raised. The name of the package from which the programs are
        available can also be included in the message.

        Internally, this program uses shutil.which, so see the documentation
        for more information about the semantics of calling.

        Arguments:
            programs (list of string): a list of programs to check

        Returns:
            list of string: a list of all programs which are not found

        Raises:
            EnvironmentError: if any programs are not callable, then
                an error is raised listing all uncallable programs.
    """
    
    import shutil

    missing_programs = []
    for program in programs:
        exe_path = shutil.which(program)

        if exe_path is None:
            missing_programs.append(program)

    if len(missing_programs) > 0:
        missing_programs = ' '.join(missing_programs)
        msg = "The following programs were not found: " + missing_programs

        if package_name is not None:
            msg = msg + ("\nPlease ensure the {} package is installed."
                .format(package_name))

        if raise_on_error:
            raise EnvironmentError(msg)
        else:
            logger.warning(msg)

    return missing_programs


def check_call_step(cmd, current_step = -1, init_step = -1, call=True, 
        raise_on_error=True):
    
    logging.info(cmd)
    ret_code = 0

    if current_step >= init_step:
        if call:
            #logging.info(cmd)
            logging.info("calling")
            ret_code = subprocess.call(cmd, shell=True)

            if raise_on_error and (ret_code != 0):
                raise subprocess.CalledProcessError(ret_code, cmd)
            elif (ret_code != 0):
                msg = ("The command returned a non-zero return code\n\t{}\n\t"
                    "Return code: {}".format(cmd, ret_code))
                logger.warning(msg)
        else:
            msg = "skipping due to --do-not-call flag"
            logging.info(msg)
    else:
        msg = "skipping due to --init-step; {}, {}".format(current_step, init_step)
        logging.info(msg)

    return ret_code


def check_call(cmd, call=True, raise_on_error=True):
    return check_call_step(cmd, call=call, raise_on_error=raise_on_error)


def check_output_step(cmd, current_step = 0, init_step = 0, raise_on_error=True):

    logging.info(cmd)
    if current_step >= init_step:
        logging.info("calling")

        try:
            out = subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError as exc:
            if raise_on_error:
                raise exc
        return out.decode()

def check_output(cmd, call=True, raise_on_error=True):
    current_step = 1
    init_step = 1
    if not call:
        init_step = 2
    return check_output_step(cmd, current_step, init_step, 
        raise_on_error=raise_on_error)


def call_if_not_exists(cmd, out_files, in_files=[], overwrite=False, call=True,
            raise_on_error=True, file_checkers=None, num_attempts=1, 
            to_delete=[], keep_delete_files=False):

    """ This function checks if out_file exists. If it does not, or if overwrite
        is true, then the command is executed, according to the call flag.
        Otherwise, a warning is issued stating that the file already exists
        and that the cmd will be skipped.

        Additionally, a list of input files can be given. If given, they must
        all exist before the call will be executed. Otherwise, a warning is 
        issued and the call is not made.

        However, if call is False, the check for input files is still made,
        but the function will continue rather than quitting. The command will
        be printed to the screen.

        The return code from the called program is returned.

        By default, if the called program returns a non-zero exit code, an
        exception is raised.

        Furthermore, a dictionary can be given which maps from a file name to
        a function which check the integrity of that file. If any of these
        function calls return False, then the relevant file(s) will be deleted
        and the call made again. The number of attempts to succeed is given as
        a parameter to the function.

        Args:
            cmd (string): the command to execute

            out_files (string or list of strings): path to the files whose existence 
                to check. If they do not exist, then the path to them will be 
                created, if necessary.

            in_files (list of strings): paths to files whose existence to check
                before executing the command

            overwrite (bool): whether to overwrite the file (i.e., execute the 
                command, even if the file exists)

            call (bool): whether to call the command, regardless of whether the
                file exists

            raise_on_error (bool): whether to raise an exception on non-zero 
                return codes

            file_checkers (dict-like): a mapping from a file name to a function
                which is used to verify that file. The function should return
                True to indicate the file is okay or False if it is corrupt. The
                functions must also accept "raise_on_error" and "logger" 
                keyword arguments.

            num_attempts (int): the number of times to attempt to create the
                output files such that all of the verifications return True.

            to_delete (list of strings): paths to files to delete if the command
                is executed successfully

            keep_delete_files (bool): if this value is True, then the to_delete
                files will not be deleted, regardless of whether the command
                succeeded

        Returns:
            int: the return code from the called program

        Warnings:
            warnings.warn if the out_file already exists and overwrite is False
            warnings.warn if the in_files do not exist

        Raises:
            subprocess.CalledProcessError: if the called program returns a
                non-zero exit code and raise_on_error is True
                
            OSError: if the maximum number of attempts is exceeded and the 
                file_checkers do not all return true and raise_on_error is True

        Imports:
            shlex
    """
    import shlex

    ret_code = 0

    # check if the input files exist
    missing_in_files = []
    for in_f in in_files:
        # we need to use shlex to ensure that we remove surrounding quotes in
        # case the file name has a space, and we are using the quotes to pass
        # it through shell
        in_f = shlex.split(in_f)[0]

        if not os.path.exists(in_f):
            missing_in_files.append(in_f)

    if len(missing_in_files) > 0:
        msg = "Some input files {} are missing. Skipping call: \n{}".format(missing_in_files, cmd)
        logger.warn(msg)
        return ret_code

        # This is here to create a directory structue using "do_not_call". In
        # hindsight, that does not seem the best way to do this, so it has been
        # removed.
        #if call:
        #    return


    # make sure we are working with a list
    if isinstance(out_files, str):
        out_files = [out_files]

    # check if the output files exist
    all_out_exists = False
    if out_files is not None:
        all_out_exists = all([os.path.exists(of) for of in out_files])

    all_valid = True
    if overwrite or not all_out_exists:
        attempt = 0
        while attempt < num_attempts:
            attempt += 1

            # create necessary paths
            if out_files is not None:
                [os.makedirs(os.path.dirname(x), exist_ok=True) for x in out_files]
            
            # make the call
            ret_code = check_call(cmd, call=call, raise_on_error=raise_on_error)

            # do not check the files if we are not calling anything
            if (not call) or (file_checkers is None):
                break

            # now check the files
            all_valid = True
            for filename, checker_function in file_checkers.items():
                msg = "Checking file for validity: {}".format(filename)
                logger.debug(msg)

                is_valid = checker_function(filename, logger=logger, 
                                raise_on_error=False)

                # if the file is not valid, then rename it
                if not is_valid:
                    all_valid = False
                    invalid_filename = "{}.invalid".format(filename)
                    msg = "Rename invalid file: {} to {}".format(filename, invalid_filename)
                    logger.warning(msg)

                    os.rename(filename, invalid_filename)

            # if they were all valid, then we are done
            if all_valid:
                break


    else:
        msg = "All output files {} already exist. Skipping call: \n{}".format(out_files, cmd)
        logger.warn(msg)

    # now, check if we succeeded in creating the output files
    if not all_valid:
        msg = ("Exceeded maximum number of attempts for cmd. The output files do "
            "not appear to be valid: {}".format(cmd))

        if raise_on_error:
            raise OSError(msg)
        else:
            logger.critical(msg)

    elif (not keep_delete_files):
        # the command succeeded, so delete the specified files
        for filename in to_delete:
            if os.path.exists(filename):
                msg = "Removing file: {}".format(filename)
                logger.info(cmd)
                
                os.remove(filename)

    return ret_code


# slurm

def check_sbatch(cmd, call=True, num_cpus=1, mem="2G", time=None, 
        partitions=None, dependencies=None, no_output=False, no_error=False, 
        use_slurm=False, mail_type=['FAIL', 'TIME_LIMIT'], mail_user=None,
        stdout_file=None, stderr_file=None,
        args=None):

    """ This function wraps calls to sbatch. It adds the relevant command line 
        options based on the parameters (either specified or extracted from 
        args, if args is not None).

        The 'ntasks' option is always 1 with the function.

        Args:

            cmd (str): The command to execute

            call (bool): If this flag is false, then the commands will not be
                executed (but will be logged).
            
            num_cpus (int): The number of CPUs to use. This will be translated into 
                an sbatch request like: "--ntasks 1 --cpus-per-task <num-cpus>". 
                default: 1

            mem (str): This will be translated into an sbatch request like: 
                "--mem=<mem>". default: 10G

            time (str): The amount of time to request. This will be translated 
                into an sbatch request like: "--time <time>". default: 0-05:59

            partitions (str): The partitions to request. This will be translated 
                into an sbatch request like: "-p <partitions>". default: general 
                (N.B. This value should be a comma-separated list with no spaces, 
                for example: partitions="general,long")

            dependencies (list of int-likes): A list of all of the job ids to
                use as dependencies for this call. This will be translated into
                an sbatch request like: "--dependency=afterok:<dependencies>".
                default: None (i.e., no dependencies)

                N.B. This IS NOT overwritten by args.
           
            no_output (bool): If this flag is True, stdout will be redirected 
                to /dev/null. This will be translated into an sbatch request 
                like: "--output=/dev/null". default: If the flag is not present, 
                then stdout will be directed to a log file with the job number. 
                This corresponds to "--output=slurm-%J.out" in the sbatch call.

            stdout_file (str): If this value is given and no_output is False,
                then this filename will be used for stdout rather than 
                slurm-%J.out. This corresponds to "--output=<stdout_file>" in 
                the sbatch call.

            no_error (bool): If this flag is True, stderr will be redirected 
                to /dev/null. This will be translated into an sbatch request 
                like: "--error=/dev/null". default: If the flag is not present, 
                then stderr will be directed to a log file with the job number. 
                This corresponds to "--error=slurm-%J.err" in the sbatch call.

            stderr_file (str): If this value is given and no_output is False,
                then this filename will be used for stderr rather than 
                slurm-%J.err. This corresponds to "--output=<stdout_file>" in 
                the sbatch call.


            use_slurm (bool): If this flag is True, then the commands will be 
                submitted to SLURM via sbatch. default: By default, each command 
                is executed sequentially within the current terminal.

            mail_type (list of strings): A list of the types of mail to send.
                This will be translated into an sbatch request like: 
                "--mail-type type_1,type_2,...". default: ['FAIL', 'TIME_LIMIT']

            mail_user (string): The email address (or user name if that is 
                configured) of the recipient of the mails. This is translated
                into an sbatch request like: "--mail-user <user>"

            args (namespace): A namespace which contains values for all of the 
                options (i.e., created from an argparse parser after calling
                add_sbatch_options on the parser)

        Returns:
            If use_slurm is False, None

            If use_slurm is True, the slurm job id
    """

    # use args if they are present
    if args is not None:
        call = not args.do_not_call
        num_cpus = args.num_cpus
        mem = args.mem
        time = args.time
        partitions = args.partitions
        no_output = args.no_output
        no_error = args.no_error
        use_slurm = args.use_slurm
        mail_type = args.mail_type
        mail_user = args.mail_user
        stdout_file = args.stdout_file
        stderr_file = args.stderr_file

    output_str = "--output=slurm-%J.out"
    if stdout_file is not None:
        output_str = "--output={}".format(stdout_file)
    if no_output:
        output_str = "--output=/dev/null"

    error_str = "--error=slurm-%J.err" 
    if stderr_file is not None:
        error_str = "--error={}".format(stderr_file)
    if no_error:
        error_str = "--error=/dev/null"

    dependencies_str = ""
    if dependencies is not None:
        dependencies_str = ':'.join(str(d) for d in dependencies)
        dependencies_str = "--dependency=afterok:{}".format(dependencies_str)

    # check if we actually want to use SLURM
    msg = "check_sbatch.use_slurm: {}, call: {}".format(use_slurm, call)
    logger.debug(msg)

    # anyway, make sure to remove the --use-slurm option
    cmd = cmd.replace("--use-slurm", "")
    
    if use_slurm:
        time_str = ""
        if time is not None:
            time_str = "--time {}".format(time)

        mem_str = ""
        if mem is not None:
            mem_str = "--mem={}".format(mem)

        partitions_str = ""
        if partitions is not None:
            partitions_str = "-p {}".format(partitions)

        num_cpus_str = ""
        if num_cpus is not None:
            num_cpus_str = "--cpus-per-task {}".format(num_cpus)

        mail_type_str = ""
        if mail_type is not None:
            mail_type_str = "--mail-type {}".format(','.join(mail_type))

        mail_user_str = ""
        if mail_user is not None:
            mail_user_str = "--mail-user {}".format(mail_user)
        else:
            # if we did not give a mail user, then do not specify the mail types
            mail_type_str = ""

        cmd = ("sbatch {} {} --ntasks 1 {}  {} "
            "{} {} {} {} {} {}".format(time_str, mem_str, partitions_str, num_cpus_str, dependencies_str, 
            output_str, error_str, mail_type_str, mail_user_str, cmd))

        output = check_output(cmd, call=call)

        # and parse out the job id
        if call:
            job_id = output.strip().split()[-1]
        else:
            job_id = None
        return job_id
    else:
        check_call(cmd, call=call)
        return None

mail_type_choices = [
    'NONE', 
    'BEGIN', 
    'END', 
    'FAIL', 
    'REQUEUE', 
    'ALL',
    'STAGE_OUT',
    'TIME_LIMIT', 
    'TIME_LIMIT_90',
    'TIME_LIMIT_80',
    'TIME_LIMIT_50',
    'ARRAY_TASKS'
]

def add_sbatch_options(parser, num_cpus=1, mem="2G", time=None, 
            stdout_file=None, stderr_file=None,
            partitions=None, mail_type=['FAIL', 'TIME_LIMIT'], mail_user=None):

    """ This function adds the options for calling sbatch to the given parser.
        The provided arguments are used as defaults for the options.

        Args:
            parser (argparse.ArgumentParser): the parser to which the
                options will be added.

            other arguments: the defaults for the sbatch options

        Returns:
            None, but parser is updated
    """
    slurm_options = parser.add_argument_group("slurm options")

    slurm_options.add_argument('--num-cpus', help="The number of CPUs to use", type=int, 
        default=num_cpus)
    slurm_options.add_argument('--mem', help="The amount of RAM to request", default=mem)
    slurm_options.add_argument('--time', help="The amount of time to request", default=time)
    slurm_options.add_argument('--partitions', help="The partitions to request", default=partitions)
    slurm_options.add_argument('--no-output', help="If this flag is present, stdout "
        "will be redirected to /dev/null", action='store_true')
    slurm_options.add_argument('--no-error', help="If this flag is present, stderr "
        "will be redirected to /dev/null", action='store_true')
    slurm_options.add_argument('--stdout-file', help="If this is present and the "
        "--no-output flag is not given, then stdout will be directed to this "
        "file rather than slurm-<job>.out", default=stdout_file)
    slurm_options.add_argument('--stderr-file', help="If this is present and the "
        "--no-error flag is not given, then stderr will be directed to this "
        "file rather than slurm-<job>.err", default=stderr_file)
    slurm_options.add_argument('--do-not-call', help="If this flag is present, then the commands "
        "will not be executed (but will be printed).", action='store_true')
    slurm_options.add_argument('--use-slurm', help="If this flag is present, the program calls "
        "will be submitted to SLURM.", action='store_true')
    slurm_options.add_argument('--mail-type', help="When to send an email notifcation "
        "of the job status. See official documentation for a description of the "
        "values. If a mail-user is not specified, this will revert to 'None'.", 
        nargs='*', choices=mail_type_choices, default=mail_type)
    slurm_options.add_argument('--mail-user', help="To whom an email will be sent, in "
        "accordance with mail-type", default=mail_user)
    
    
# utils


def abspath(*fn):
    return os.path.abspath(os.path.join(os.sep, *fn))


def check_keys_exist(d, keys):
    """ This function ensures the given keys are present in the dictionary. It
        does not other validate the type, value, etc., of the keys or their
        values. If a key is not present, a KeyError is raised.

        The motivation behind this function is to verify that a config dictionary
        read in at the beginning of a program contains all of the required values.
        Thus, the program will immediately detect when a required config value is
        not present and quit.

        Input:
            d (dict) : the dictionary

            keys (list) : a list of keys to check
        Returns:
            list of string: a list of all programs which are not found

        Raises:
            KeyError: if any of the keys are not in the dictionary
    """
    missing_keys = [k for k in keys if k not in d]

    
    if len(missing_keys) > 0:
        missing_keys = ' '.join(missing_keys)
        msg = "The following keys were not found: " + missing_keys
        raise KeyError(msg)

    return missing_keys


def check_files_exist(files, raise_on_error=True, logger=logger, 
        msg="The following files were missing: ", source=None):
    """ This function ensures that all of the files in the list exists. If any
        do not, it will either raise an exception or print a warning, depending
        on the value of raise_on_error.

        Parameters
        ----------
        files: list of strings
            the file paths to check

        raise_on_error: bool
            whether to raise an error if any of the files are missing

        logger: logging.Logger
            a logger to use for writing the warning if an error is not raised

        msg: string
            a message to write before the list of missing files

        source: string
            a description of where the check is made, such as a module name. If
            this is not None, it will be prepended in brackets before msg.

        Returns
        -------
        all_exist: bool
            True if all of the files existed, False otherwise

        Raises
        ------
        FileNotFoundError, if raise_on_error is True and any of the files
                do not exist.
    """
    missing_files = []

    for f in files:
        if not os.path.exists(f):
            missing_files.append(f)

    if len(missing_files) == 0:
        return True

    missing_files_str = ",".join(missing_files)
    source_str = ""
    if source is not None:
        source_str = "[{}]: ".format(source)
    msg = "{}{}{}".format(source_str, msg, missing_files_str)
    if raise_on_error:
        raise FileNotFoundError(msg)
    else:
        logger.warning(msg)

    return False


# parallel

def apply_parallel_iter(items, num_procs, func, *args, progress_bar=False, total=None, num_groups=None, backend='loky'):
    """ This function parallelizes applying a function to all items in an iterator using the 
        joblib library. In particular, func is called for each of the items in the list. (Unless
        num_groups is given. In this case, the iterator must be "split-able", e.g., a list or
        np.array. Then func is called with a list of items.)
        
        This function is best used when func has little overhead compared to the item processing.

        It always returns a list of the values returned by func. The order of the
        returned list is dependent on the semantics of joblib.Parallel, but it is typically
        in the same order as the groups

        Args:
            items (list-like): A list (or anything that can be iterated over in a for loop)

            num_procs (int): The number of processors to use

            func (function pointer): The function to apply to each group

            args (variable number of arguments): The other arguments to pass to func

            progress_bar (boolean) : whether to use a tqdm progress bar

            total (int) : the total number of items in the list. This only needs to be used if
                len(items) is not defined. It is only used for tqdm, so it is not necessary.

            num_groups (int) : if given, the number of groups into which the input is split. In
                case it is given, then each call to func will be passed a list of items.

        Returns:
            list: the values returned from func for each item (in the order specified by
                joblib.Parallel). If num_groups is given, then this is likely to be a list.

        Imports:
            joblib
            numpy
            tqdm, if progress_bar is True
    """
    import joblib
    import numpy as np

    # check if we want groups
    if num_groups is not None:
        items_per_group = int(np.ceil(len(items) / num_groups))
        # then make items a list of lists, where each internal list contains items from the original list
        items = [ items[i * items_per_group : (i+1) * items_per_group] for i in range(num_groups)]

    if progress_bar:
        import tqdm
        ret_list = joblib.Parallel(n_jobs=num_procs, backend=backend)(joblib.delayed(func)(item, *args) 
            for item in tqdm.tqdm(items, leave=True, file=sys.stdout, total=total))
    else:
        ret_list = joblib.Parallel(n_jobs=num_procs, backend=backend)(joblib.delayed(func)(item, *args) for item in items)
    return ret_list


# bam_utils


def check_bam_file(filename, check_index=False, raise_on_error=True, logger=logger):
    """ This function wraps a call to "samtools view" and pipes the output to 
        /dev/null. Additionally, it optionally checks that the index is present.
        
        Optionally, it raises an exception if the return code is 
        not 0. Otherwise, it writes a "critical" warning message.

        Args:
            filename (str): a path to the bam file

            check_index (bool): check whether the index is present

            raise_on_error (bool): whether to raise an OSError (if True) or log
                a "critical" message (if false)

            logger (logging.Logger): a logger for writing the message if an
                error is not raised

        Returns:
            bool: whether the file was valid

        Raises:
            OSError: if quickcheck does not return 0 and raise_on_error is True

    """
    
    programs = ['samtools']
    check_programs_exist(programs)

    dev_null = abspath("dev", "null")

    cmd = "samtools view -h {} > {}".format(filename, dev_null)

    ret = check_call_step(cmd, raise_on_error=False)

    if ret != 0:
        msg = "The bam/sam file does not appear to be valid: {}".format(filename)

        if raise_on_error:
            raise OSError(msg)

        logger.critical(msg)
        return False

    # now look for the index
    if not check_index:
        return True

    cmd = "samtools idxstats {}".format(filename)

    ret = check_call_step(cmd, raise_on_error=False)

    if ret != 0:
        msg = "The bam/sam file does not appear to have a valid index: {}".format(filename)

        if raise_on_error:
            raise OSError(msg)

        logger.critical(msg)
        return False


    # then the file and the index was okay
    return True


def sort_bam_file(bam, sorted_bam, args):
    """ Sort bam file, wrapping a call to "samtools sort" via
    call_if_not_exists.

    Args:
        bam (str): path to bam file
        sorted_bam (str): path to sorted bam file
        args (namespace): calling arguments
    """

    call = not args.do_not_call
    keep_delete_files = args.keep_intermediate_files or args.do_not_call

    sam_tmp_str = ""
    if args.tmp is not None:
        sam_tmp_dir = os.path.join(args.tmp, "{}_samtools".format(args.name))
        if not os.path.exists(sam_tmp_dir):
            os.makedirs(sam_tmp_dir)
        sam_tmp_str = "-T {}".format(sam_tmp_dir)

    cmd = "samtools sort {} -@{} -o {} {}".format(
        bam,
        args.num_cpus,
        sorted_bam,
        sam_tmp_str
    )

    in_files = [bam]
    out_files = [sorted_bam]
    to_delete = [bam]
    file_checkers = {
        sorted_bam: check_bam_file
    }
    call_if_not_exists(cmd, out_files, in_files=in_files,
        file_checkers=file_checkers, overwrite=args.overwrite, call=call,
        keep_delete_files=keep_delete_files, to_delete=to_delete)


def index_bam_file(bam, args):
    """ Index bam file, wrapping a call to "samtools index" via
    call_if_not_exists.

    Args:
        bam (str): path to bam file
        args (namespace): calling arguments
    """

    call = not args.do_not_call

    bai = bam + ".bai"
    cmd = "samtools index -b {} {}".format(bam, bai)
    in_files = [bam]
    out_files = [bai]
    call_if_not_exists(cmd, out_files, in_files=in_files,
        overwrite=args.overwrite, call=call)
    
    
def get_pysam_alignment_file(f, mode=None, **kwargs):
    """ This function checks the type of a given object and returns a pysam
        AlignmentFile object. If the object is already an AlignmentFile, it
        is simply returned. If it is a string, then it is treated as a file
        path and opened as an AlignmentFile. Otherwise, an Exception is raised.

        Args:
            f (obj): either an existing AlignmentFile or a path to a bam/sam file

            mode (str): the mode for opening the file, if f is a file path

            **kwargs : other keyword arguments to pass to the AlignmentFile
                constructor

        Returns:
            pysam.AlignmentFile: the AlignmentFile for f

        Raises:
            TypeError: if f is neither an AlignmentFile nor string
        Import:
            pysam
    """

    import pysam

    if isinstance(f, pysam.AlignmentFile):
        return f

    elif isinstance(f, str):
        return pysam.AlignmentFile(f, mode=mode, **kwargs)

    else:
        msg = "Could not interpret value as pysam.AlignmentFile: {}".format(f)
        raise ValueError(msg)
    
    

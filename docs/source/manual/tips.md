# Tips & Tricks

## Useful Snakemake CLI options

```
usage: snakemake [OPTIONS] -- [TARGET ...]

  --dry-run, --dryrun, -n
                        Do not execute anything, and display what would be done. If you have a very large workflow,
                        use `--dry-run --quiet` to just print a summary of the DAG of jobs. (default: False)
  --workflow-profile WORKFLOW_PROFILE
                        Path (relative to current directory) to workflow specific profile folder to use for
                        configuring Snakemake with parameters specific for this workflow (like resources).
  --jobs, -j N          Use at most N CPU cluster/cloud jobs in parallel. For local execution this is an alias for
                        `--cores` (it is though recommended to use `--cores` in that case). Note: Set to `unlimited`
                        to allow any number of parallel jobs.
  --config, -C [KEY=VALUE ...]
                        Set or overwrite values in the workflow config object. The workflow config object is
                        accessible as variable config inside the workflow. Default values can be set by providing a
                        YAML JSON file (see `--configfile` and Documentation).
  --keep-going, -k      Go on with independent jobs if a job fails during execution. This only applies to runtime
                        failures in job execution, not to errors during workflow parsing or DAG construction.
                        (default: False)
  --rerun-triggers {code,input,mtime,params,software-env} [{code,input,mtime,params,software-env} ...]
                        Define what triggers the rerunning of a job. By default, all triggers are used, which
                        guarantees that results are consistent with the workflow code and configuration. If you
                        rather prefer the traditional way of just considering file modification dates, use `--rerun-
                        trigger mtime`. (default: code input mtime params software-env)
  --rerun-triggers {code,input,mtime,params,software-env} [{code,input,mtime,params,software-env} ...]
                        Define what triggers the rerunning of a job. By default, all triggers are used, which
                        guarantees that results are consistent with the workflow code and configuration. If you
                        rather prefer the traditional way of just considering file modification dates, use `--rerun-
                        trigger mtime`. (default: code input mtime params software-env)
  --force, -f           Force the execution of the selected target or the first rule regardless of already created
                        output. (default: False)
  --forceall, -F        Force the execution of the selected (or the first) rule and all rules it is dependent on
                        regardless of already created output. (default: False)
  --forcerun, -R [TARGET ...]
                        Force the re-execution or creation of the given rules or files. Use this option if you
                        changed a rule and want to have all its output in your workflow updated.
  --rerun-incomplete, --ri
                        Re-run all jobs the output of which is recognized as incomplete. (default: False)
  --list-rules, --list, -l
                        Show available rules in given Snakefile. (default: False)
  --summary, -S         Print a summary of all files created by the workflow. The has the following columns:
                        filename, modification time, rule version, status, plan. Thereby rule version contains the
                        version the file was created with (see the version keyword of rules), and status denotes
                        whether the file is missing, its input files are newer or if version or implementation of
                        the rule changed since file creation. Finally the last column denotes whether the file will
                        be updated or created during the next workflow execution. (default: False)
```

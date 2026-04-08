def write2bdfoutput(bdfname, grpid, data_ses):
    jk = len(data_ses[0])  # Get current length of data

    key1 = f'jobfile.open( "{bdfname}", "ANALYZE NO JOBFILE" )'
    key2 = f'msc_delete_old_files( "{bdfname}", ".bdf", ".op2" )'
    key3 = f'jobfile.writec( "JOBNAME", "{bdfname}" )'
    key4 = f'jobfile.writec( "GROUP", "{grpid}" )'
    key5 = 'jobfile.writei( "SELECTED GROUP 0", 1 )'
    key6 = f'jobfile.writec( "SELECTED GROUP 1", "{grpid}" )'
    key7 = f'mscnastran_job.associate_subcases( "101", "{bdfname}", 1, ["Default"] )'
    key8 = f'analysis_submit_2( "MSC.Nastran", "{bdfname}" )'

    # Read session file
    with open("group2bdf.ses", "r") as f:
        data3 = f.readlines()
    # Inject keys at specific lines, rest are original
    for I, line in enumerate(data3):
        if I == 0:
            data_ses[0].append(key1)
        elif I == 1:
            data_ses[0].append(key2)
        elif I == 7:
            data_ses[0].append(key3)
        elif I == 18:
            data_ses[0].append(key4)
        elif I == 19:
            data_ses[0].append(key5)
        elif I == 20:
            data_ses[0].append(key6)
        elif I == 170:
            data_ses[0].append(key7)
        elif I == 171:
            data_ses[0].append(key8)
        else:
            data_ses[0].append(line.strip())  # Strip newline for consistency

    return data_ses

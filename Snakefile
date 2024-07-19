rule key:
    output:
        "psg_status.txt"
    script:
        "src/scripts/setup_psg.py"

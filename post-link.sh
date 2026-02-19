#!/bin/bash

# ==============================================================================
# GMAAMP Post-Installation Script
# ==============================================================================

LOG_FILE="$PREFIX/.messages.txt"
echo "=== GenoMycAnalyzer Post-Installation Setup ===" >> "$LOG_FILE"

# ------------------------------------------------------------------------------
# 1. SnpEff DB Configuration
# ------------------------------------------------------------------------------
GMA_DB_DIR="${PREFIX}/share/GMA/snpEff_DB/data/H37Rv"
SNPEFF_JAR=$(find ${PREFIX}/share -name "snpEff.jar" | head -n 1)

if [ -n "$SNPEFF_JAR" ]; then
    SNPEFF_HOME=$(dirname "$SNPEFF_JAR")
    mkdir -p "${SNPEFF_HOME}/data"
    
    if [ -d "${SNPEFF_HOME}/data/H37Rv" ] && [ ! -L "${SNPEFF_HOME}/data/H37Rv" ]; then
        mv "${SNPEFF_HOME}/data/H37Rv" "${SNPEFF_HOME}/data/H37Rv_backup_pkg"
    fi
    
    rm -f "${SNPEFF_HOME}/data/H37Rv"
    
    ln -sfn "$GMA_DB_DIR" "${SNPEFF_HOME}/data/H37Rv"
    
    SNPEFF_CONFIG="${SNPEFF_HOME}/snpEff.config"
    if ! grep -q "Custom MTB H37Rv" "$SNPEFF_CONFIG"; then
        cat <<EOF >> "$SNPEFF_CONFIG"

# --- GMA Custom MTB H37Rv ---
H37Rv.genome : Mycobacterium tuberculosis H37Rv
	H37Rv.chromosomes : NC_000962.3
	H37Rv.NC_000962.3.codonTable : Bacterial_and_Plant_Plastid
EOF
    fi
    echo " [Success] SnpEff database configured." >> "$LOG_FILE"
else
    echo " [Warning] SnpEff.jar not found. Skipping configuration." >> "$LOG_FILE"
fi

# ------------------------------------------------------------------------------
# 2. TB-Profiler Database Update
# ------------------------------------------------------------------------------
TBP_EXE="$PREFIX/bin/tb-profiler"
if [ -x "$TBP_EXE" ]; then
    echo " [Info] Updating TB-Profiler database... (This may take a moment)" >> "$LOG_FILE"
    "$TBP_EXE" update_tbdb >> "$LOG_FILE" 2>&1 || true
    echo " [Success] TB-Profiler database update step finished." >> "$LOG_FILE"
fi

# ------------------------------------------------------------------------------
# 3. Kraken2 Database Notice
# ------------------------------------------------------------------------------
echo "" >> "$LOG_FILE"
echo " [IMPORTANT] Kraken2 Database Requirement" >> "$LOG_FILE"
echo " GMA requires a Kraken2 database to run. This package does NOT include it" >> "$LOG_FILE"
echo " due to its large size. Please prepare a database using:" >> "$LOG_FILE"
echo "   $ kraken2-build --standard --threads 16 --db YOUR_DB_NAME" >> "$LOG_FILE"
echo " or download a pre-built index." >> "$LOG_FILE"
echo "" >> "$LOG_FILE"

echo "=== Setup Complete ===" >> "$LOG_FILE"
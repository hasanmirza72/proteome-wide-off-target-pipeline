import os
import glob
import re

def main():
    print("====================================================================")
    print("ALL SCRIPTS in scripts/  (name, size, and files each references)")
    print("====================================================================")

    # Get all .py and .sh scripts in the scripts/ folder
    scripts = sorted(glob.glob("scripts/*.py") + glob.glob("scripts/*.sh"))
    
    # Regex to find data file names enclosed in quotes
    file_pattern = re.compile(r'["\']([A-Za-z0-9_./-]+\.(?:tsv|csv|json|txt|fasta|pdb))["\']')

    for script in scripts:
        try:
            with open(script, 'r', encoding='utf-8', errors='ignore') as f:
                lines = f.readlines()
        except Exception:
            continue
            
        line_count = len(lines)
        print(f"\n### {script}  ({line_count} lines)")

        # Find purpose (the very first line starting with # or """)
        purpose = ""
        for line in lines[:20]:
            stripped = line.strip()
            if stripped.startswith('#') or stripped.startswith('"""') or stripped.startswith("'''"):
                purpose = stripped
                break
        
        if purpose:
            print(f"   purpose: {purpose}")
        else:
            print("   purpose: (none found)")

        # Find referenced files (limit to top 12 alphabetically, just like the bash script)
        content = "".join(lines)
        references = set(file_pattern.findall(content))
        
        print("   references:")
        if references:
            for ref in sorted(references)[:12]:
                print(f"      {ref}")
        else:
            print("      (none)")

    print("\n====================================================================")
    print(f"COUNT: {len(scripts)} scripts total in scripts/")
    print("====================================================================")
    
    print("\nAlso list any scripts OUTSIDE scripts/ (root, tools/, extra_files/):")
    
    # Check for stray scripts in other folders
    outside_patterns = [
        "*.py", "*.sh", 
        "tools/*.py", "tools/*.sh", 
        "extra_files/*.py", "extra_files/*.sh"
    ]
    
    outside_scripts = []
    for pattern in outside_patterns:
        outside_scripts.extend(glob.glob(pattern))
        
    # Remove any scripts from the scripts/ folder just in case
    outside_scripts = [s for s in outside_scripts if not s.startswith("scripts/")]
    
    if outside_scripts:
        for s in sorted(outside_scripts):
            print(f"  - {s}")
    else:
        print("  (none in root/tools/extra_files)")

if __name__ == "__main__":
    main()

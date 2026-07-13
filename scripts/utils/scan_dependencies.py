import os
import ast
import re

# Standard Python libraries we want to ignore so they don't clutter your list
STD_LIB = {
    'os', 'sys', 'csv', 're', 'json', 'time', 'datetime', 'math', 'shutil', 
    'subprocess', 'argparse', 'collections', 'itertools', 'typing', 'urllib',
    'multiprocessing', 'concurrent', 'glob', 'io', 'random', 'logging', 'requests'
}

def main():
    scripts_dir = 'scripts'
    if not os.path.exists(scripts_dir):
        print(f"Error: Could not find the folder '{scripts_dir}'")
        return

    python_modules = set()
    urls = set()
    bash_commands = set()

    # Regex to find links (http/https)
    url_pattern = re.compile(r'https?://[^\s\'"<>]+')

    for root, _, files in os.walk(scripts_dir):
        for file in files:
            path = os.path.join(root, file)
            
            # --- 1. SCAN PYTHON FILES ---
            if file.endswith('.py'):
                with open(path, 'r', encoding='utf-8', errors='ignore') as f:
                    content = f.read()
                    
                    # Grab URLs
                    urls.update(url_pattern.findall(content))
                    
                    # Parse Python imports safely
                    try:
                        tree = ast.parse(content)
                        for node in ast.walk(tree):
                            if isinstance(node, ast.Import):
                                for alias in node.names:
                                    base_mod = alias.name.split('.')[0]
                                    if base_mod not in STD_LIB:
                                        python_modules.add(base_mod)
                            elif isinstance(node, ast.ImportFrom):
                                if node.module:
                                    base_mod = node.module.split('.')[0]
                                    if base_mod not in STD_LIB:
                                        python_modules.add(base_mod)
                    except SyntaxError:
                        print(f"Skipping syntax error in {file}")

            # --- 2. SCAN BASH FILES ---
            elif file.endswith('.sh'):
                with open(path, 'r', encoding='utf-8', errors='ignore') as f:
                    content = f.read()
                    
                    # Grab URLs
                    urls.update(url_pattern.findall(content))
                    
                    # Guess command line tools (first word of a line)
                    for line in content.split('\n'):
                        line = line.strip()
                        if line and not line.startswith('#'):
                            # Ignore basic bash keywords
                            first_word = line.split()[0]
                            if first_word not in ('echo', 'cd', 'export', 'if', 'for', 'while', 'done', 'fi', 'then', 'else', 'rm', 'mkdir', 'ls', 'cat', 'cp', 'mv'):
                                bash_commands.add(first_word)

    # --- PRINT THE RESULTS ---
    print("\n========== DEPENDENCY REPORT ==========")
    
    print("\n📦 EXTERNAL PYTHON MODULES USED:")
    if python_modules:
        for mod in sorted(python_modules):
            print(f"  - {mod}")
    else:
        print("  (None found)")

    print("\n🔗 URLs AND LINKS FOUND:")
    if urls:
        for url in sorted(urls):
            print(f"  - {url}")
    else:
        print("  (None found)")

    print("\n🛠️  EXTERNAL CLI TOOLS CALLED (from .sh files):")
    if bash_commands:
        for cmd in sorted(bash_commands):
            print(f"  - {cmd}")
    else:
        print("  (None found)")
    print("\n=======================================\n")

if __name__ == "__main__":
    main()

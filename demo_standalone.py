"""
Quantum Drug Discovery - Standalone Demo
Run complete pipeline locally with toy Hamiltonian, print ΔG, save plots
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), 'qdrug', 'src'))

import numpy as np

try:
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    plt = None  # type: ignore
    rcParams = None  # type: ignore
    print("Warning: matplotlib not available. Install with `pip install matplotlib` to enable plotting.")

# Configure plotting when matplotlib is available
if MATPLOTLIB_AVAILABLE:
    rcParams['figure.dpi'] = 100
    rcParams['font.size'] = 11
    rcParams['axes.grid'] = True
    rcParams['grid.alpha'] = 0.3
else:
    print("Running in text-only mode; plots will be skipped.")

from parse_pdb import create_h2_molecule, create_lih_molecule, FragmentExtractor
from hamiltonian import HamiltonianBuilder
from vqe_runner import VQERunner
from binding_calc import BindingCalculator
from utils import hartree_to_kj
import io
import html
import webbrowser

def print_header(title):
    """Print formatted section header"""
    print("\n" + "="*70)
    print(f"  {title}")
    print("="*70)

# Capture terminal output
class OutputCapture:
    def __init__(self):
        self.output = []
    
    def write(self, text):
        self.output.append(text)
        sys.__stdout__.write(text)
    
    def flush(self):
        sys.__stdout__.flush()
    
    def get_output(self):
        return ''.join(self.output)

def plot_pes_curve(bond_lengths, vqe_energies, exact_energies, filename='pes_curve.png'):
    """Plot potential energy surface"""
    if not MATPLOTLIB_AVAILABLE:
        print("  ⚠️ matplotlib not installed; skipping PES plot (run `pip install matplotlib`).")
        return
    plt.figure(figsize=(10, 6))
    plt.plot(bond_lengths, exact_energies, 'o-', label='Exact', 
             linewidth=2.5, markersize=8, color='#004AAD')
    plt.plot(bond_lengths, vqe_energies, 's--', label='VQE', 
             linewidth=2, markersize=6, color='#0073E5', alpha=0.8)
    plt.xlabel('Bond Length (Å)', fontsize=13, fontweight='bold')
    plt.ylabel('Energy (Hartree)', fontsize=13, fontweight='bold')
    plt.title('H₂ Potential Energy Surface', fontsize=15, fontweight='bold')
    plt.legend(fontsize=12, framealpha=0.9)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {filename}")
    plt.close()

def plot_vqe_convergence(iterations, energies, final_energy, filename='vqe_convergence.png'):
    """Plot VQE convergence"""
    if not MATPLOTLIB_AVAILABLE:
        print("  ⚠️ matplotlib not installed; skipping convergence plot (run `pip install matplotlib`).")
        return
    plt.figure(figsize=(10, 6))
    plt.plot(iterations, energies, 'o-', linewidth=2, markersize=4, color='#0073E5')
    plt.axhline(y=final_energy, color='#FF4444', linestyle='--', 
                linewidth=2, label=f"Final: {final_energy:.6f} Ha")
    plt.xlabel('Iteration', fontsize=13, fontweight='bold')
    plt.ylabel('Energy (Hartree)', fontsize=13, fontweight='bold')
    plt.title('VQE Convergence for H₂', fontsize=15, fontweight='bold')
    plt.legend(fontsize=12, framealpha=0.9)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {filename}")
    plt.close()

def plot_binding_comparison(results, filename='binding_comparison.png'):
    """Plot comprehensive ΔG bar chart"""
    if not MATPLOTLIB_AVAILABLE:
        print("  ⚠️ matplotlib not installed; skipping binding comparison plot (run `pip install matplotlib`).")
        return
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
    
    names = list(results.keys())
    delta_e_elec = [results[n]['delta_e_electronic_kj_mol'] for n in names]
    delta_g_solv = [results[n]['delta_g_solvation_kj_mol'] for n in names]
    t_delta_s = [results[n]['t_delta_s_kj_mol'] for n in names]
    
    x = np.arange(len(names))
    width = 0.65
    
    # Stacked components
    ax1.bar(x, delta_e_elec, width, label='ΔE_elec', color='#0073E5')
    ax1.bar(x, delta_g_solv, width, bottom=delta_e_elec, 
            label='ΔG_solv', color='#4DA8FF')
    ax1.bar(x, t_delta_s, width, 
            bottom=np.array(delta_e_elec)+np.array(delta_g_solv), 
            label='-TΔS', color='#CCE5FF')
    
    ax1.set_xlabel('System', fontsize=13, fontweight='bold')
    ax1.set_ylabel('Energy (kJ/mol)', fontsize=13, fontweight='bold')
    ax1.set_title('ΔG Components Breakdown', fontsize=15, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(names, rotation=15, ha='right')
    ax1.legend(fontsize=11, framealpha=0.9)
    ax1.axhline(y=0, color='black', linestyle='-', linewidth=1)
    ax1.grid(alpha=0.3)
    
    # Total ΔG with color coding
    delta_g_total = [results[n]['delta_g_binding_kj_mol'] for n in names]
    colors = ['#00AA00' if dg < -30 else '#FFA500' if dg < -20 else '#FF4444' 
              for dg in delta_g_total]
    
    bars = ax2.bar(x, delta_g_total, width, color=colors, 
                   edgecolor='black', linewidth=1.5, alpha=0.9)
    ax2.axhline(y=-30, color='green', linestyle='--', linewidth=2, 
                label='Strong (<-30 kJ/mol)', alpha=0.7)
    ax2.axhline(y=-20, color='orange', linestyle='--', linewidth=2, 
                label='Moderate (<-20 kJ/mol)', alpha=0.7)
    
    ax2.set_xlabel('System', fontsize=13, fontweight='bold')
    ax2.set_ylabel('ΔG_binding (kJ/mol)', fontsize=13, fontweight='bold')
    ax2.set_title('Total Binding Free Energy', fontsize=15, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(names, rotation=15, ha='right')
    ax2.legend(fontsize=11, framealpha=0.9)
    ax2.grid(alpha=0.3)
    
    # Value labels
    for bar, val in zip(bars, delta_g_total):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                 f'{val:.1f}', ha='center', 
                 va='bottom' if val < 0 else 'top', 
                 fontweight='bold', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {filename}")
    plt.close()

def generate_html_report(terminal_output, best_system, best_dg, h2_error):
    """Generate HTML report with terminal logs and visualizations"""
    
    # Escape and format terminal output for HTML
    escaped_output = html.escape(terminal_output)
    
    # Add color classes to terminal output
    escaped_output = escaped_output.replace('Warning:', '<span class="warning">Warning:</span>')
    escaped_output = escaped_output.replace('✓', '<span class="success">✓</span>')
    escaped_output = escaped_output.replace('✗', '<span class="error">✗</span>')
    escaped_output = escaped_output.replace('converged SCF energy', '<span class="value">converged SCF energy</span>')
    
    # Highlight section headers
    lines = escaped_output.split('\n')
    formatted_lines = []
    for line in lines:
        if '=====' in line or 'PART' in line or 'QUANTUM DRUG DISCOVERY' in line or 'SUMMARY' in line:
            formatted_lines.append(f'<span class="header-line">{line}</span>')
        elif 'Drug-' in line or 'Reference:' in line:
            formatted_lines.append(f'<span class="header-line">{line}</span>')
        elif 'Ha' in line or 'kJ/mol' in line or 'nM' in line or 'Qubits:' in line or 'Electrons:' in line:
            # Highlight numerical values
            import re
            line = re.sub(r'(-?\d+\.?\d*(?:e[+-]?\d+)?)\s*(Ha|kJ/mol|nM)', 
                         r'<span class="value">\1 \2</span>', line)
            line = re.sub(r'(Qubits|Electrons|Iterations|Time):\s*(\d+)', 
                         r'\1: <span class="value">\2</span>', line)
            formatted_lines.append(line)
        else:
            formatted_lines.append(line)
    
    escaped_output = '\n'.join(formatted_lines)
    
    html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Quantum Drug Discovery - Results</title>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        :root {{
            --primary: #667eea;
            --primary-dark: #5568d3;
            --secondary: #764ba2;
            --success: #10b981;
            --warning: #f59e0b;
            --error: #ef4444;
            --dark: #1f2937;
            --light: #f9fafb;
            --border: #e5e7eb;
            --shadow: rgba(0, 0, 0, 0.1);
        }}
        
        body {{
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            padding: 20px;
            color: var(--dark);
        }}
        
        .container {{
            max-width: 1600px;
            margin: 0 auto;
            background: white;
            border-radius: 24px;
            box-shadow: 0 25px 50px -12px rgba(0, 0, 0, 0.25);
            overflow: hidden;
        }}
        
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 60px 40px;
            text-align: center;
            position: relative;
            overflow: hidden;
        }}
        
        .header::before {{
            content: '';
            position: absolute;
            top: -50%;
            left: -50%;
            width: 200%;
            height: 200%;
            background: radial-gradient(circle, rgba(255,255,255,0.1) 1px, transparent 1px);
            background-size: 50px 50px;
            animation: drift 20s linear infinite;
            pointer-events: none;
        }}
        
        @keyframes drift {{
            0% {{ transform: translate(0, 0); }}
            100% {{ transform: translate(50px, 50px); }}
        }}
        
        .header-content {{
            position: relative;
            z-index: 1;
        }}
        
        .header h1 {{
            font-size: 3em;
            margin-bottom: 15px;
            font-weight: 800;
            letter-spacing: -0.02em;
            text-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        }}
        
        .header p {{
            font-size: 1.3em;
            opacity: 0.95;
            font-weight: 400;
        }}
        
        .nav-tabs {{
            display: flex;
            background: #f8f9fa;
            border-bottom: 2px solid var(--border);
            position: sticky;
            top: 0;
            z-index: 100;
        }}
        
        .nav-tab {{
            flex: 1;
            padding: 24px;
            text-align: center;
            cursor: pointer;
            font-weight: 600;
            color: #6b7280;
            transition: all 0.3s ease;
            border: none;
            background: transparent;
            font-size: 1.05em;
            position: relative;
        }}
        
        .nav-tab::after {{
            content: '';
            position: absolute;
            bottom: -2px;
            left: 0;
            right: 0;
            height: 3px;
            background: var(--primary);
            transform: scaleX(0);
            transition: transform 0.3s ease;
        }}
        
        .nav-tab:hover {{
            background: #e5e7eb;
            color: var(--primary);
        }}
        
        .nav-tab.active {{
            background: white;
            color: var(--primary);
        }}
        
        .nav-tab.active::after {{
            transform: scaleX(1);
        }}
        
        .content {{
            padding: 50px;
        }}
        
        .tab-content {{
            display: none;
            animation: fadeInUp 0.6s ease;
        }}
        
        .tab-content.active {{
            display: block;
        }}
        
        @keyframes fadeInUp {{
            from {{ 
                opacity: 0; 
                transform: translateY(20px);
            }}
            to {{ 
                opacity: 1; 
                transform: translateY(0);
            }}
        }}
        
        .terminal {{
            background: #0d1117;
            color: #c9d1d9;
            padding: 40px;
            border-radius: 16px;
            font-family: 'JetBrains Mono', 'Courier New', monospace;
            font-size: 13px;
            line-height: 1.8;
            overflow-x: auto;
            white-space: pre-wrap;
            word-wrap: break-word;
            max-height: 700px;
            overflow-y: auto;
            box-shadow: inset 0 2px 10px rgba(0, 0, 0, 0.5), 0 4px 20px rgba(0, 0, 0, 0.3);
            border: 1px solid #30363d;
        }}
        
        .terminal::-webkit-scrollbar {{
            width: 12px;
        }}
        
        .terminal::-webkit-scrollbar-track {{
            background: #161b22;
            border-radius: 8px;
        }}
        
        .terminal::-webkit-scrollbar-thumb {{
            background: #30363d;
            border-radius: 8px;
            border: 2px solid #161b22;
        }}
        
        .terminal::-webkit-scrollbar-thumb:hover {{
            background: #484f58;
        }}
        
        .terminal .success {{
            color: #3fb950;
            font-weight: 500;
        }}
        
        .terminal .warning {{
            color: #d29922;
        }}
        
        .terminal .error {{
            color: #f85149;
            font-weight: 500;
        }}
        
        .terminal .header-line {{
            color: #58a6ff;
            font-weight: 600;
        }}
        
        .terminal .value {{
            color: #7ee787;
        }}
        
        .visualization {{
            margin-bottom: 60px;
            background: white;
            border-radius: 16px;
            overflow: hidden;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08);
            border: 1px solid var(--border);
            transition: all 0.3s ease;
        }}
        
        .visualization:hover {{
            box-shadow: 0 8px 30px rgba(0, 0, 0, 0.12);
            transform: translateY(-2px);
        }}
        
        .visualization h2 {{
            color: var(--dark);
            padding: 30px 30px 20px;
            font-size: 1.8em;
            font-weight: 700;
            background: linear-gradient(135deg, #f8f9fa 0%, #e5e7eb 100%);
            border-bottom: 3px solid var(--primary);
        }}
        
        .visualization img {{
            width: 100%;
            display: block;
            transition: transform 0.3s ease;
        }}
        
        .visualization:hover img {{
            transform: scale(1.01);
        }}
        
        .viz-description {{
            padding: 30px;
            color: #6b7280;
            line-height: 1.8;
            background: #fafafa;
            font-size: 1.05em;
        }}
        
        .summary-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            border-radius: 20px;
            margin-bottom: 40px;
            box-shadow: 0 10px 40px rgba(102, 126, 234, 0.3);
            position: relative;
            overflow: hidden;
        }}
        
        .summary-card::before {{
            content: '';
            position: absolute;
            top: 0;
            right: 0;
            width: 300px;
            height: 300px;
            background: radial-gradient(circle, rgba(255,255,255,0.1) 0%, transparent 70%);
            border-radius: 50%;
            transform: translate(30%, -30%);
        }}
        
        .summary-card h3 {{
            font-size: 1.8em;
            margin-bottom: 30px;
            font-weight: 700;
            position: relative;
        }}
        
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
            gap: 25px;
            position: relative;
        }}
        
        .summary-item {{
            background: rgba(255, 255, 255, 0.15);
            padding: 30px;
            border-radius: 16px;
            backdrop-filter: blur(20px);
            border: 1px solid rgba(255, 255, 255, 0.2);
            transition: all 0.3s ease;
        }}
        
        .summary-item:hover {{
            background: rgba(255, 255, 255, 0.2);
            transform: translateY(-4px);
            box-shadow: 0 8px 20px rgba(0, 0, 0, 0.15);
        }}
        
        .summary-item .label {{
            font-size: 0.95em;
            opacity: 0.9;
            margin-bottom: 12px;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            font-weight: 500;
        }}
        
        .summary-item .value {{
            font-size: 2em;
            font-weight: 800;
            letter-spacing: -0.02em;
        }}
        
        .section-title {{
            font-size: 2em;
            font-weight: 700;
            margin-bottom: 30px;
            color: var(--dark);
            display: flex;
            align-items: center;
            gap: 12px;
        }}
        
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 25px;
            margin-top: 30px;
        }}
        
        .stat-card {{
            background: white;
            padding: 35px 30px;
            border-radius: 16px;
            text-align: center;
            border: 2px solid var(--border);
            transition: all 0.3s ease;
            position: relative;
            overflow: hidden;
        }}
        
        .stat-card::before {{
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            height: 4px;
            background: linear-gradient(90deg, var(--primary), var(--secondary));
            transform: scaleX(0);
            transition: transform 0.3s ease;
        }}
        
        .stat-card:hover {{
            border-color: var(--primary);
            transform: translateY(-4px);
            box-shadow: 0 8px 30px rgba(102, 126, 234, 0.15);
        }}
        
        .stat-card:hover::before {{
            transform: scaleX(1);
        }}
        
        .stat-card .stat-value {{
            font-size: 2.5em;
            font-weight: 800;
            background: linear-gradient(135deg, var(--primary), var(--secondary));
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
            margin-bottom: 12px;
            letter-spacing: -0.02em;
        }}
        
        .stat-card .stat-label {{
            color: #6b7280;
            font-size: 1em;
            font-weight: 500;
        }}
        
        .download-section {{
            margin-top: 50px;
            padding: 40px;
            background: linear-gradient(135deg, #f8f9fa 0%, #e5e7eb 100%);
            border-radius: 20px;
            border: 2px dashed var(--border);
        }}
        
        .download-btn {{
            display: inline-block;
            background: linear-gradient(135deg, var(--primary), var(--secondary));
            color: white;
            padding: 16px 36px;
            border-radius: 12px;
            text-decoration: none;
            font-weight: 600;
            font-size: 1.05em;
            transition: all 0.3s ease;
            margin: 10px 12px 10px 0;
            box-shadow: 0 4px 15px rgba(102, 126, 234, 0.3);
            border: none;
            cursor: pointer;
        }}
        
        .download-btn:hover {{
            transform: translateY(-2px);
            box-shadow: 0 8px 25px rgba(102, 126, 234, 0.4);
        }}
        
        .download-btn:active {{
            transform: translateY(0);
        }}
        
        .refresh-btn {{
            background: linear-gradient(135deg, #10b981, #059669);
        }}
        
        .refresh-btn:hover {{
            box-shadow: 0 8px 25px rgba(16, 185, 129, 0.4);
        }}
        
        @media (max-width: 768px) {{
            .header h1 {{
                font-size: 2em;
            }}
            
            .content {{
                padding: 30px 20px;
            }}
            
            .nav-tab {{
                padding: 16px;
                font-size: 0.9em;
            }}
            
            .stats-grid, .summary-grid {{
                grid-template-columns: 1fr;
            }}
        }}
        
        .badge {{
            display: inline-block;
            padding: 6px 14px;
            border-radius: 20px;
            font-size: 0.85em;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.05em;
        }}
        
        .badge-success {{
            background: #d1fae5;
            color: #065f46;
        }}
        
        .badge-info {{
            background: #dbeafe;
            color: #1e40af;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="header-content">
                <h1>⚛️ Quantum Drug Discovery</h1>
                <p>Comprehensive Analysis & Visualization Results</p>
            </div>
        </div>
        
        <div class="nav-tabs">
            <button class="nav-tab active" onclick="showTab('summary')">📊 Summary</button>
            <button class="nav-tab" onclick="showTab('visualizations')">📈 Visualizations</button>
            <button class="nav-tab" onclick="showTab('logs')">🖥️ Terminal Logs</button>
        </div>
        
        <div class="content">
            <div id="summary" class="tab-content active">
                <div class="summary-card">
                    <h3>🎯 Execution Summary</h3>
                    <div class="summary-grid">
                        <div class="summary-item">
                            <div class="label">Status</div>
                            <div class="value">✓ Complete</div>
                        </div>
                        <div class="summary-item">
                            <div class="label">Best Binder</div>
                            <div class="value">{best_system}</div>
                        </div>
                        <div class="summary-item">
                            <div class="label">ΔG Binding</div>
                            <div class="value">{best_dg:.2f} kJ/mol</div>
                        </div>
                        <div class="summary-item">
                            <div class="label">Systems Evaluated</div>
                            <div class="value">5</div>
                        </div>
                    </div>
                </div>
                
                <h2 class="section-title">📊 Key Metrics</h2>
                <div class="stats-grid">
                    <div class="stat-card">
                        <div class="stat-value">H₂</div>
                        <div class="stat-label">Molecule Validated</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-value">4</div>
                        <div class="stat-label">Qubits Used</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-value">8</div>
                        <div class="stat-label">Bond Lengths Scanned</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-value">{h2_error:.6f}</div>
                        <div class="stat-label">VQE Error (Ha)</div>
                    </div>
                </div>
                
                <div class="download-section">
                    <h2 class="section-title">📥 Download Results</h2>
                    <a href="pes_curve.png" download class="download-btn">⬇️ PES Curve</a>
                    <a href="vqe_convergence.png" download class="download-btn">⬇️ VQE Convergence</a>
                    <a href="binding_comparison.png" download class="download-btn">⬇️ Binding Analysis</a>
                    <button onclick="location.reload()" class="download-btn refresh-btn">🔄 Refresh</button>
                </div>
            </div>
            
            <div id="visualizations" class="tab-content">
                <div class="visualization">
                    <h2>1️⃣ Potential Energy Surface (PES)</h2>
                    <img src="pes_curve.png" alt="PES Curve">
                    <div class="viz-description">
                        <strong>Analysis:</strong> The potential energy surface demonstrates how molecular energy varies with bond length. 
                        The comparison between exact (classical) and VQE (quantum) calculations validates the accuracy 
                        of the quantum algorithm across different molecular configurations.
                    </div>
                </div>
                
                <div class="visualization">
                    <h2>2️⃣ VQE Convergence Analysis</h2>
                    <img src="vqe_convergence.png" alt="VQE Convergence">
                    <div class="viz-description">
                        <strong>Analysis:</strong> This convergence plot tracks the Variational Quantum Eigensolver (VQE) optimization process 
                        over multiple iterations, illustrating how the algorithm efficiently converges to the ground state 
                        energy of the H₂ molecule.
                    </div>
                </div>
                
                <div class="visualization">
                    <h2>3️⃣ Binding Free Energy Comparison</h2>
                    <img src="binding_comparison.png" alt="Binding Comparison">
                    <div class="viz-description">
                        <strong>Analysis:</strong> Comprehensive binding energy analysis featuring component breakdown (electronic energy, 
                        solvation effects, entropy) and total ΔG values for multiple drug candidates. More negative ΔG values 
                        indicate stronger binding affinity and better potential therapeutic efficacy.
                    </div>
                </div>
            </div>
            
            <div id="logs" class="tab-content">
                <h2 class="section-title">🖥️ Terminal Output</h2>
                <div class="terminal">
{escaped_output}
                </div>
            </div>
        </div>
    </div>
    
    <script>
        function showTab(tabName) {{
            const tabs = document.querySelectorAll('.tab-content');
            tabs.forEach(tab => tab.classList.remove('active'));
            
            const navTabs = document.querySelectorAll('.nav-tab');
            navTabs.forEach(tab => tab.classList.remove('active'));
            
            document.getElementById(tabName).classList.add('active');
            event.target.classList.add('active');
        }}
    </script>
</body>
</html>'''
    
    return html_content

def main():
    """Run complete quantum drug discovery demo"""
    
    # Capture output
    output_capture = OutputCapture()
    original_stdout = sys.stdout
    sys.stdout = output_capture
    
    print("\n" + "="*70)
    print("  QUANTUM DRUG DISCOVERY - STANDALONE DEMO")
    print("="*70)
    print("\n  This demo runs with toy/mock Hamiltonians (no PySCF/Qiskit required)")
    print("  Generates: PES curve, VQE convergence, ΔG comparison charts\n")
    
    # Part 1: H2 Validation
    print_header("PART 1: H₂ MOLECULE VALIDATION")
    
    h2_geom = create_h2_molecule(bond_length=0.74)
    builder = HamiltonianBuilder(basis='sto-3g')
    h2_ham = builder.build_from_geometry(h2_geom)
    
    print(f"\nMolecule: H₂")
    print(f"  Qubits: {h2_ham['num_qubits']}")
    print(f"  Electrons: {h2_ham['num_electrons']}")
    print(f"  HF Energy: {h2_ham['hf_energy']:.6f} Ha")
    
    runner = VQERunner()
    comparison = runner.compare_with_exact(h2_ham)
    
    print(f"\nVQE Validation:")
    print(f"  Exact: {comparison['exact_energy']:.6f} Ha")
    print(f"  VQE:   {comparison['vqe_energy']:.6f} Ha")
    print(f"  Error: {comparison['error']:.6f} Ha ({comparison['error_percent']:.2f}%)")
    print(f"  Chemical Accuracy: {'✓ YES' if comparison['within_chemical_accuracy'] else '✗ NO'}")
    
    # Part 2: Potential Energy Surface
    print_header("PART 2: POTENTIAL ENERGY SURFACE")
    
    bond_lengths = np.linspace(0.5, 2.0, 8)
    vqe_energies = []
    exact_energies = []
    
    print("\nScanning H₂ bond lengths...")
    for r in bond_lengths:
        geom = create_h2_molecule(bond_length=r)
        ham = builder.build_from_geometry(geom)
        comp = runner.compare_with_exact(ham)
        vqe_energies.append(comp['vqe_energy'])
        exact_energies.append(comp['exact_energy'])
        print(f"  r = {r:.2f} Å: VQE = {comp['vqe_energy']:.6f} Ha")
    
    plot_pes_curve(bond_lengths, vqe_energies, exact_energies)
    
    # Part 3: VQE Convergence
    print_header("PART 3: VQE CONVERGENCE TRACKING")
    
    h2_geom = create_h2_molecule()
    h2_ham = builder.build_from_geometry(h2_geom)
    
    print("\nRunning VQE with convergence tracking...")
    result = runner.run_vqe(h2_ham, reps=2, max_iterations=40)
    
    history = result['convergence_history']
    iterations = [h['iteration'] for h in history]
    energies = [h['energy'] for h in history]
    
    print(f"  Final Energy: {result['energy']:.6f} Ha")
    print(f"  Iterations: {result['optimizer_evals']}")
    print(f"  Time: {result['optimizer_time']:.2f} s")
    
    plot_vqe_convergence(iterations, energies, result['energy'])
    
    # Part 4: Binding Energy Calculations
    print_header("PART 4: BINDING ENERGY CALCULATIONS")
    
    systems = {
        'Drug-A': {'complex': -22.5, 'protein': -16.0, 'ligand': -5.8},
        'Drug-B': {'complex': -20.2, 'protein': -15.0, 'ligand': -4.9},
        'Drug-C': {'complex': -24.1, 'protein': -17.2, 'ligand': -6.2},
        'Drug-D': {'complex': -19.5, 'protein': -14.5, 'ligand': -4.7},
        'Reference': {'complex': -18.0, 'protein': -13.8, 'ligand': -4.0},
    }
    
    calculator = BindingCalculator(temperature=298.15)
    results = {}
    
    print("\nCalculating ΔG for multiple systems:")
    for name, energies in systems.items():
        result = calculator.calculate_binding_energy(
            {'energy': energies['complex']},
            {'energy': energies['protein']},
            {'energy': energies['ligand']},
            include_corrections=True
        )
        results[name] = result
        
        print(f"\n{name}:")
        print(f"  E_complex:  {energies['complex']:.2f} Ha")
        print(f"  E_protein:  {energies['protein']:.2f} Ha")
        print(f"  E_ligand:   {energies['ligand']:.2f} Ha")
        print(f"  ΔE_elec:    {result['delta_e_electronic_kj_mol']:+.2f} kJ/mol")
        print(f"  ΔG_solv:    {result['delta_g_solvation_kj_mol']:+.2f} kJ/mol")
        print(f"  -TΔS:       {result['t_delta_s_kj_mol']:+.2f} kJ/mol")
        print(f"  ──────────────────────────")
        print(f"  ΔG_binding: {result['delta_g_binding_kj_mol']:+.2f} kJ/mol")
        
        # IC50 estimate
        ic50 = calculator.calculate_ic50_estimate(result['delta_g_binding_kj_mol'])
        print(f"  IC₅₀ est:   {ic50:.2e} nM")
    
    plot_binding_comparison(results)
    
    # Summary
    print_header("SUMMARY")
    
    delta_g_values = [results[n]['delta_g_binding_kj_mol'] for n in results.keys()]
    best_system = list(results.keys())[np.argmin(delta_g_values)]
    best_dg = min(delta_g_values)
    
    print(f"\n✓ Complete pipeline executed successfully!")
    print(f"\nGenerated Files:")
    if MATPLOTLIB_AVAILABLE:
        print(f"  • pes_curve.png (Potential Energy Surface)")
        print(f"  • vqe_convergence.png (VQE Optimization)")
        print(f"  • binding_comparison.png (ΔG Bar Charts)")
    else:
        print("  • Plot generation skipped (install matplotlib to enable saved figures)")
    print(f"\nKey Results:")
    print(f"  • H₂ VQE Error: {comparison['error']:.6f} Ha ({comparison['error_percent']:.2f}%)")
    print(f"  • Best Binder: {best_system} (ΔG = {best_dg:.2f} kJ/mol)")
    print(f"  • Systems Evaluated: {len(systems)}")
    
    print("\n" + "="*70)
    print("  Demo complete! Check generated PNG files for visualizations.")
    print("="*70 + "\n")
    
    # Restore stdout
    sys.stdout = original_stdout
    
    # Generate HTML report
    terminal_output = output_capture.get_output()
    html_report = generate_html_report(terminal_output, best_system, best_dg, comparison['error'])
    
    # Save HTML report
    html_filename = 'results.html'
    with open(html_filename, 'w', encoding='utf-8') as f:
        f.write(html_report)
    
    print(f"✓ Generated HTML report: {html_filename}")
    print(f"  Opening in browser...")
    
    # Open in default browser
    import os
    webbrowser.open('file://' + os.path.abspath(html_filename))

if __name__ == "__main__":
    main()

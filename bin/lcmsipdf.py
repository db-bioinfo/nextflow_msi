import json
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.platypus import Image, PageBreak
import argparse
from datetime import datetime

def create_title_style():
    styles = getSampleStyleSheet()
    title_style = ParagraphStyle(
        'CustomTitle',
        parent=styles['Heading1'],
        fontSize=16,
        spaceAfter=30,
        alignment=1  # Center alignment
    )
    return title_style

def create_header_style():
    styles = getSampleStyleSheet()
    header_style = ParagraphStyle(
        'CustomHeader',
        parent=styles['Heading2'],
        fontSize=12,
        spaceAfter=12,
        textColor=colors.HexColor('#2E5C8A')
    )
    return header_style

def create_normal_style():
    styles = getSampleStyleSheet()
    normal_style = ParagraphStyle(
        'CustomNormal',
        parent=styles['Normal'],
        fontSize=10,
        spaceAfter=6
    )
    return normal_style

def create_status_table(data):
    status_color = colors.red if data['msi_status']['classification'] == 'MSI-H' else colors.green
    status_data = [
        ['MSI Status', data['msi_status']['classification']],
        ['Unstable Sites', f"{data['msi_status']['unstable_sites']}/{data['msi_status']['total_sites']} ({data['msi_status']['instability_percentage']:.2f}%)"],
        ['Quality Status', data['quality_status']['overall']],
        ['Clinical Implication', 'Potential response to immunotherapy' if data['msi_status']['classification'] == 'MSI-H' else 'Standard treatment protocols recommended']
    ]
    
    table = Table(status_data, colWidths=[2*inch, 4*inch])
    table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (0, -1), colors.HexColor('#F5F5F5')),
        ('TEXTCOLOR', (1, 0), (1, 0), status_color),
        ('TEXTCOLOR', (1, 2), (1, 2), 
         colors.green if data['quality_status']['overall'] == 'PASS' else colors.red),
        ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
        ('FONTNAME', (0, 0), (-1, -1), 'Helvetica'),
        ('FONTSIZE', (0, 0), (-1, -1), 10),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
    ]))
    return table

def create_quality_metrics_table(data):
    metrics = data['quality_metrics']['all_sites']
    min_coverage = data['technical_details']['coverage_thresholds']['minimum']
    optimal_coverage = data['technical_details']['coverage_thresholds']['optimal']
    
    qc_data = [
        ['Quality Metric', 'Value', 'Threshold', 'Status'],
        ['Mean Coverage', f"{metrics['mean_coverage']:.1f}x", f"≥{min_coverage}x", 
         'PASS' if metrics['mean_coverage'] >= min_coverage else 'FAIL'],
        ['Median Coverage', f"{metrics['median_coverage']:.1f}x", f"≥{min_coverage}x",
         'PASS' if metrics['median_coverage'] >= min_coverage else 'FAIL'],
        ['Sites Above Min Coverage', f"{(metrics['sites_above_min_coverage']/metrics['total_sites']*100):.1f}%", "≥95%",
         'PASS' if metrics['sites_above_min_coverage']/metrics['total_sites'] >= 0.95 else 'FAIL'],
        ['Sites Above Optimal', f"{(metrics['sites_above_optimal']/metrics['total_sites']*100):.1f}%", f"≥{optimal_coverage}x",
         'PASS' if metrics['sites_above_optimal']/metrics['total_sites'] >= 0.80 else 'WARNING']
    ]
    
    table = Table(qc_data, colWidths=[2*inch, 1.5*inch, 1.5*inch, inch])
    table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#F5F5F5')),
        ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
        ('FONTNAME', (0, 0), (-1, -1), 'Helvetica'),
        ('FONTSIZE', (0, 0), (-1, -1), 10),
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
        ('TEXTCOLOR', (3, 1), (3, -1), 
         colors.green if data['quality_status']['overall'] == 'PASS' else colors.red),
    ]))
    return table

def create_notable_sites_table(data):
    if not data['notable_sites']:
        return None
        
    table_data = [['Location', 'Repeat Unit', 'Times', 'Instability', 'Coverage', 'Status']]
    for site in data['notable_sites']:
        status_color = colors.green if site['quality_status'] == 'PASS' else colors.red
        table_data.append([
            site['location'],
            site['repeat_unit'],
            str(site['repeat_times']),
            f"{site['instability']*100:.1f}%",
            f"{site['coverage']}x",
            site['quality_status']
        ])
    
    table = Table(table_data, colWidths=[1.5*inch, inch, 0.7*inch, inch, inch, 0.8*inch])
    table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#F5F5F5')),
        ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
        ('FONTNAME', (0, 0), (-1, -1), 'Helvetica'),
        ('FONTSIZE', (0, 0), (-1, -1), 9),
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
    ]))
    return table

def create_warnings_section(data, normal_style):
    if not data['quality_status'].get('warnings') and not data['quality_status'].get('flags'):
        return None
    
    elements = []
    if data['quality_status'].get('flags'):
        elements.append(Paragraph("Quality Flags:", create_header_style()))
        for flag in data['quality_status']['flags']:
            elements.append(Paragraph(f"• {flag}", normal_style))
        elements.append(Spacer(1, 12))
    
    if data['quality_status'].get('warnings'):
        elements.append(Paragraph("Warnings:", create_header_style()))
        for warning in data['quality_status']['warnings']:
            elements.append(Paragraph(f"• {warning}", normal_style))
    
    return elements

def generate_report(json_file, output_pdf):
    # Load the analysis data
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Create the document
    doc = SimpleDocTemplate(
        output_pdf,
        pagesize=letter,
        rightMargin=72,
        leftMargin=72,
        topMargin=72,
        bottomMargin=72
    )
    
    # Styles
    title_style = create_title_style()
    header_style = create_header_style()
    normal_style = create_normal_style()
    
    # Build the document
    elements = []
    
    # Title and Header
    elements.append(Paragraph(f"MSI Analysis Report", title_style))
    elements.append(Paragraph(f"Sample ID: {data['sample_id']}", header_style))
    elements.append(Paragraph(f"Analysis Date: {data['analysis_date']}", normal_style))
    elements.append(Spacer(1, 20))
    
    # MSI Status Section
    elements.append(Paragraph("MSI Status and Quality", header_style))
    elements.append(create_status_table(data))
    elements.append(Spacer(1, 20))
    
    # Quality Metrics Section
    elements.append(Paragraph("Quality Metrics", header_style))
    elements.append(create_quality_metrics_table(data))
    elements.append(Spacer(1, 20))
    
    # Warnings Section (if any)
    warnings_elements = create_warnings_section(data, normal_style)
    if warnings_elements:
        elements.extend(warnings_elements)
        elements.append(Spacer(1, 20))
    
    # Notable Sites Section
    if data['notable_sites']:
        elements.append(Paragraph("Notable Unstable Sites", header_style))
        elements.append(create_notable_sites_table(data))
        elements.append(Spacer(1, 20))
    
    # Technical Details
    elements.append(Paragraph("Technical Details", header_style))
    tech_details = f"""
    Analysis performed using {data['technical_details']['tool']}
    Minimum coverage threshold: {data['technical_details']['coverage_thresholds']['minimum']}x
    Optimal coverage threshold: {data['technical_details']['coverage_thresholds']['optimal']}x
    MSI threshold: {data['technical_details']['instability_threshold']*100}%
    """
    elements.append(Paragraph(tech_details, normal_style))
    
    # Generate the PDF
    doc.build(elements)

def main():
    parser = argparse.ArgumentParser(description='Generate MSI Analysis PDF Report')
    parser.add_argument('--input', required=True, help='Input JSON analysis file')
    parser.add_argument('--output', required=True, help='Output PDF file path')
    
    args = parser.parse_args()
    
    try:
        generate_report(args.input, args.output)
        print(f"Report generated: {args.output}")
    except Exception as e:
        print(f"Error generating report: {str(e)}")
        raise

if __name__ == "__main__":
    main()

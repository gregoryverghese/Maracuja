
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Maracuja README</title>
</head>
<body>
    <h1>Maracuja <img src="maracuja.png" alt="Maracuja Logo"></h1>
    <p>Maracuja is a comprehensive package for parsing and analyzing TCR sequencing data. It supports data from various platforms, providing tools for detailed clonality analysis.</p>

    <h2>Table of Contents</h2>
    <ul>
        <li><a href="#parsing">Parsing</a></li>
        <li><a href="#analysis">Analysis</a>
            <ul>
                <li><a href="#clonality">Clonality</a>
                    <ul>
                        <li><a href="#stacked-bar-chart">Stacked Bar Chart</a></li>
                        <li><a href="#grouped-bar-chart">Grouped Bar Chart</a></li>
                    </ul>
                </li>
            </ul>
        </li>
    </ul>

    <h2 id="parsing">Parsing</h2>
    <p>Maracuja supports parsing raw TCR sequencing files from Adaptive Biotechnologies and ImmunoSEQ.</p>

    <h3>Supported Formats</h3>
    <ul>
        <li>Adaptive Biotechnologies</li>
        <li>ImmunoSEQ</li>
    </ul>

    <h3>Parsing Example</h3>
    <p>To parse the raw files, use the <code>rep_load</code> function:</p>
    <pre><code>folder_path = 'path/to/your/folder_with_repertoire_files/'
data, meta_df = rep_load(folder_path)</code></pre>

    <h2 id="analysis">Analysis</h2>
    <p>Maracuja provides powerful tools to analyze clonality within TCR sequencing data.</p>

    <h3 id="clonality">Clonality</h3>

    <h4 id="stacked-bar-chart">Stacked Bar Chart</h4>
    <p>You can create a stacked bar chart to visualize the top proportions of clones within your samples.</p>
    <pre><code>top_prop = top_proportion(data)
plot_top_proportion_stacked(top_prop)</code></pre>
    <img src="top.png" alt="Stacked Bar Chart">

    <h4 id="grouped-bar-chart">Grouped Bar Chart</h4>
    <p>Stratify your samples based on a particular group to create a grouped bar chart using rare proportions.</p>
    <pre><code>rare_prop = rare_proportion(data)
plot_top_proportion_grouped(rare_prop, meta_df, 'Group_Column_Name')</code></pre>
    <img src="rare.png" alt="Grouped Bar Chart">

    <h2>Installation</h2>
    <p>To install Maracuja, use:</p>
    <pre><code>pip install maracuja</code></pre>

    <h2>Usage</h2>
    <p>After installation, import the package and use the provided functions to parse and analyze your TCR sequencing data.</p>

    <h2>Contributing</h2>
    <p>We welcome contributions! Please fork the repository and submit a pull request.</p>

    <h2>License</h2>
    <p>Maracuja is licensed under the MIT License.</p>
</body>
</html>


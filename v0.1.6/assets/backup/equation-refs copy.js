// Automatically number equation references based on CSS counter
document.addEventListener('DOMContentLoaded', function() {
    // Wait for KaTeX rendering and CSS counters to be applied
    setTimeout(function() {
        // Build a map of equation IDs to their numbers
        const equationMap = {};
        const h6Elements = document.querySelectorAll('.docs-main h6');
        
        h6Elements.forEach((h6, index) => {
            const anchor = h6.querySelector('a.docs-heading-anchor');
            if (anchor) {
                const id = anchor.getAttribute('href').substring(1); // Remove leading #
                equationMap[id] = index + 1; // CSS counter starts at 1
            }
        });
        
        // Update all links that reference equations (but NOT the ones in h6 headers)
        const links = document.querySelectorAll('a[href^="#"]');
        links.forEach(link => {
            // Skip if this link is inside an h6 header
            if (link.closest('h6')) {
                return;
            }
            
            const href = link.getAttribute('href').substring(1);
            if (equationMap[href] !== undefined) {
                const text = link.textContent.trim();
                // Check if it's an equation reference pattern
                if (text.startsWith('Eq.') || text.match(/^Eq\.\s*\d+/)) {
                    // Extract the label name if present
                    const match = text.match(/Eq\.\s*\d*:?\s*(.+)/);
                    if (match && match[1]) {
                        link.textContent = `Eq. ${equationMap[href]}: ${match[1]}`;
                    } else {
                        link.textContent = `Eq. ${equationMap[href]}`;
                    }
                } else if (!text.match(/^\d+$/)) {
                    // If it's just the label name, prepend equation number
                    link.textContent = `Eq. ${equationMap[href]}: ${text}`;
                }
            }
        });
    }, 500); // Delay to ensure CSS counters are applied
});


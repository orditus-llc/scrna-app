// Placeholder drawing of the Orditus mark (concentric green rings with a light spark).
// Replace with the official logo file when available (CLAUDE.md, "UI design").
export function Logo({ size = 44 }: { size?: number }) {
  return (
    <svg width={size} height={size} viewBox="0 0 44 44" role="img" aria-label="Orditus logo">
      <circle cx="22" cy="22" r="21" fill="#CFE8DA" />
      <circle cx="22" cy="22" r="17" fill="#8FC2A9" />
      <circle cx="22" cy="22" r="13" fill="#4E8571" />
      <circle cx="22" cy="22" r="9" fill="#1F4A3F" />
      <circle cx="22" cy="22" r="5.5" fill="#0B1F1C" />
      <path d="M22 8 L22.9 21.5 L22 23.5 L21.1 21.5 Z" fill="#F4FBF7" />
      <circle cx="22" cy="22.6" r="1.5" fill="#FFFFFF" />
    </svg>
  );
}

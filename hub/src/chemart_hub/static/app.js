// Chemart Hub: the only script. Copy buttons; everything else works without it.
document.addEventListener("click", async (event) => {
  const button = event.target.closest("[data-copy]");
  if (!button) return;
  const source = document.getElementById(button.dataset.copy);
  const text = source ? (source.value ?? source.innerText) : button.dataset.copy;
  try {
    await navigator.clipboard.writeText(text.trim());
    const label = button.textContent;
    button.textContent = "Copied";
    setTimeout(() => { button.textContent = label; }, 1400);
  } catch (err) {
    if (source && source.select) source.select();
  }
});

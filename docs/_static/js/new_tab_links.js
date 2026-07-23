document.addEventListener("DOMContentLoaded", () => {
  document.querySelectorAll(".new-tab a").forEach((link) => {
    link.target = "_blank";
    link.rel = "noopener";
  });
});
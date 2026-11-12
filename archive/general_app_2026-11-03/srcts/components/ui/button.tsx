// shadcn/ui-style button, restyled with the Option B tokens.
import { cva, type VariantProps } from "class-variance-authority";
import type { ButtonHTMLAttributes } from "react";

import { cn } from "../../lib/utils";

const buttonVariants = cva(
  "inline-flex items-center justify-center gap-2 rounded-[4px] text-sm font-bold transition-colors " +
    "focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-ink " +
    "disabled:cursor-not-allowed disabled:opacity-100",
  {
    variants: {
      variant: {
        primary: "bg-ink text-white hover:bg-[#2a2a2e] disabled:bg-line disabled:text-muted",
        secondary:
          "border border-field bg-white text-ink enabled:hover:bg-chip disabled:border-line disabled:bg-chip disabled:text-off",
      },
      size: {
        default: "px-[18px] py-[10px]",
        small: "px-3 py-[7px] text-[13px]",
      },
    },
    defaultVariants: { variant: "primary", size: "default" },
  },
);

export interface ButtonProps
  extends ButtonHTMLAttributes<HTMLButtonElement>,
    VariantProps<typeof buttonVariants> {}

export function Button({ className, variant, size, type = "button", ...props }: ButtonProps) {
  return <button type={type} className={cn(buttonVariants({ variant, size }), className)} {...props} />;
}

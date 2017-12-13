function print_figure(h,output_dir,filename)

print(h,[output_dir filename '.png'],'-dpng','-r600');
print(h,[output_dir filename '.eps'],'-depsc');

end
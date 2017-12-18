require 'rubygems'
require 'time'
require 'sequel'
#require 'database.rb'
require 'ogs_author_mapping.rb'

Sequel::Model.unrestrict_primary_key

# Create a table if it not exists
$DB.create_table? :commit_infos do
  primary_key     :revision
  Time            :date
  String          :branch
  foreign_key     :author_id, :table => :authors
end

class CommitInfo < Sequel::Model(:commit_infos)
  set_primary_key :revision
  set_dataset dataset.order(:revision)
  many_to_one :author
  many_to_one :benchmark_run
end

class CommitInfoLoader

  def new?
    @new
  end

  def initialize(filename)

    @new = true

    File.open(filename, 'r') do |file|
      revision = 0
      author = nil
      date = nil
      branch = nil

      while line = file.gets
        line.scan(/svn\/ogs\/([\S]+)\//) do |match|
          branch = match[0]
        end
        line.scan(/Revision:\s([0-9]+)/) do |match|
          revision = match[0].to_i
        end
        line.scan(/Last Changed Author:\s([\w]+)/) do |match|
          author_name = match[0]
          author = Author[:svn_user => author_name]
        end
        line.scan(/Last Changed Date:\s([0-9]{4}-[0-9]{2}-[0-9]{2}\s[0-9]{2}:[0-9]{2}:[0-9]{2})/) do |match|
          date = Time.parse(match[0])
        end
      end

      if CommitInfo[:revision => revision]
        @new = false
        puts "Commit info of revision #{revision} already read."
      else
        commit_info = CommitInfo.create(:revision => revision,
                                        :date => date,
                                        :branch => branch)
        commit_info.author = author
        commit_info.save
      end

    end
  end

end

#CommitInfoLoader.new('tests/svnInfoOld.txt')
#$DB[:commit_infos].each {|row| p row}